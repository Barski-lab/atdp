#include "atdpbasics.hpp"


ATDPBasics::ATDPBasics(EXPERIMENT_INFO* e)
{
    exp_i=e;
    avd_window=gArgs().getArgs("avd_window").toInt();
    avd_whole_region=avd_window*2+1;
    avd_heat_window=gArgs().getArgs("avd_heat_window").toInt();
    twicechr=gArgs().getArgs("sam_twicechr").toString().toLower();
    ignorechr=gArgs().getArgs("sam_ignorechr").toString().toLower();
    avd_bodysize=1000;

    /* Open bam files preparing EXPERIMENT_INFO class
     */
    QString fp = exp_i->filepath;
    if ( !exp_i->reader.Open(fp.toStdString()) ) {
        throw "Could not open input BAM files";
    }

    exp_i->header = exp_i->reader.GetHeader();
    exp_i->references = exp_i->reader.GetReferenceData();

    //Organazing chrom ids
    for(int RefID=0;RefID < (int)exp_i->references.size();RefID++) {
        if(ignorechr.contains(QString(exp_i->references[RefID].RefName.c_str()).toLower())) {
            exp_i->i_tids<<RefID;
            continue;
        }
        if(twicechr.contains(QString(exp_i->references[RefID].RefName.c_str()).toLower())) {
            exp_i->tids<<RefID;
        }
        exp_i->ref_map.insert(
                    exp_i->references[RefID].RefName.c_str(),
                    QPair<int,int>(RefID,exp_i->references[RefID].RefLength));
    }

    if(gArgs().getArgs("index").toString().isEmpty()) {
        if (not exp_i->reader.LocateIndex(BamIndex::STANDARD)){
            throw "Couldn't locate index alongside the BAM file";
        };
    } else if (not exp_i->reader.OpenIndex(gArgs().getArgs("index").toString().toStdString())) {
        throw "Could not locate index set by --index parameter";
    };

    //Fill in all regions for the current uid
    this->getRegions();
}

void ATDPBasics::getRegions() {

    if(gArgs().getArgs("annotation").toString().isEmpty()) {
        throw "Set gene annotaion file";
    }

    string annotation_filename = gArgs().getArgs("annotation").toString().toStdString();
    ifstream annotation_stream (annotation_filename.c_str());
    if (!annotation_stream) {
        throw "Failed to open annotaion file";
    }

    std::string line;
    getline(annotation_stream, line);
    std::map<std::string,short> column_map;
    if (line.length() > 0 && line.at(0) == '#'){
        QStringList header_splitted = QString::fromStdString(line).split(QChar('\t'));
        for (int i = 0; i < header_splitted.length(); i++){
            column_map[header_splitted.at(i).toStdString()] = i;
        }
    } else {
        column_map["exonCount"] = 8;
        column_map["exonStarts"] = 9;
        column_map["exonEnds"] = 10;
        column_map["name"] = 1;
        column_map["name2"] = 12;
        column_map["chrom"] = 2;
        column_map["strand"] = 3;
        column_map["txStart"] = 4;
        column_map["txEnd"] = 5;
    }

    while(getline(annotation_stream, line)) {
        QStringList line_splitted = QString::fromStdString(line).split(QChar('\t'));

        QString chr = line_splitted.at(column_map["chrom"]);
        QChar strand = line_splitted.at(column_map["strand"]).at(0);
        quint64 txStart = line_splitted.at(column_map["txStart"]).toInt(); // Do we need to add +1
        quint64 txEnd = line_splitted.at(column_map["txEnd"]).toInt();
        QString refseq_id = line_splitted.at(column_map["name"]);
        QString gene_id = line_splitted.at(column_map["name2"]);

        if(ignorechr.contains(chr)) {
            continue;
        }

        QSharedPointer<REGION> region(new REGION());
        qint64 start=0,end=0;

        if(strand=='+'){
            start=txStart-avd_window;//-exp_i->fragmentsize/2;
            end=txStart+avd_window;
            region->strand=true;
        } else {
            start=txEnd-avd_window;//-exp_i->fragmentsize/2;
            end=txEnd+avd_window;
            region->strand=false;
        }

        if(!exp_i->regions.isEmpty() && exp_i->regions.last()->start==start) {
            if(!exp_i->regions.last()->gene_id.contains(gene_id))
                exp_i->regions.last()->gene_id.append(","+gene_id);
            if(!exp_i->regions.last()->refseq_id.contains(refseq_id))
                exp_i->regions.last()->refseq_id.append(","+refseq_id);
            region.clear();
            continue;
        }

        region->txStart=txStart;
        region->txEnd=txEnd;

        region->start=start;
        region->end=end;//+exp_i->fragmentsize/2;

        region->gene_id=gene_id;
        region->refseq_id=refseq_id;
        region->chrom=chr;
        exp_i->avd_matrix.append(
                    QPair<QSharedPointer<REGION>,QVector<quint16> >(
                        QSharedPointer<REGION>(region),QVector<quint16>(avd_whole_region/avd_heat_window,0) ) );
        QJsonArray body;
        for(int c=0;c<300;c++)
            body<<0.0;
        exp_i->body_matrix.append(
                    QPair<QSharedPointer<REGION>,QJsonArray >(
                        QSharedPointer<REGION>(region), body) );

        exp_i->regions.append(region);
    }
}

void ATDPBasics::RegionsProcessing () {
    /*
     *  Work with all regions
     */
    double r_w_b=(double)avd_window/(double)avd_bodysize;

    for(int i=0;i<exp_i->regions.size();i++) {
        if(!exp_i->reader.SetRegion(
               exp_i->ref_map[exp_i->regions[i]->chrom].first,
               exp_i->regions[i]->txStart-avd_window,
               exp_i->ref_map[exp_i->regions[i]->chrom].first,
               exp_i->regions[i]->txEnd+avd_window
               )){
            qDebug()<<"Cant set region:"
                   <<exp_i->regions[i]->chrom
                  <<exp_i->regions[i]->start
                 <<exp_i->regions[i]->end;
            continue;
            //throw "Cant set region.";
        } else {
            BamAlignment al;
            int shift=exp_i->fragmentsize/2;
            qint64 left=exp_i->regions[i]->start;//+shift;

            while ( exp_i->reader.GetNextAlignmentCore(al) ) {

                if(al.IsPaired() && (!al.IsProperPair())) {
                    continue;
                }

                if(al.IsPaired() && al.IsProperPair() && al.IsMateMapped() ) {
                    if( ( (al.Position<al.MatePosition) && al.IsReverseStrand() ) || ( (al.MatePosition < al.Position) && al.IsMateReverseStrand() )) {
                        prn_debug("Name2:",al);
                        continue;
                    }
                }

                int position_b = al.Position+1;
                int position_e = al.GetEndPosition();

                int length=abs(al.InsertSize);

                if(al.IsMateMapped() && al.IsFirstMate() ) { //pair-end reads
                    if( length==0 ) { //bug
                        prn_debug("Name3:",al);
                        continue;
                    }
                    if(al.IsReverseStrand()) {
                        position_b= position_e-length+1;
                    } else {
                        position_e= position_b+length-1;
                    }
                } else  if(al.IsMateMapped() && al.IsSecondMate()) {
                    continue;
                }

                int b;//Most left position for Average tag density
                int b1;//Current read position shifted to fragment size for Average Gene Body
                if(al.IsPaired()) {
                    b=position_b-left+length/2;
                    b1=position_b+length/2;
                } else {
                    if(al.IsReverseStrand()) {// - strand
                        b=position_e-shift-left;
                        b1=position_e-shift;
                    } else {
                        b=position_b+shift-left;
                        b1=position_b+shift;
                    }
                }

                // Doubles requested chromosomes
                quint16 d=1;
                if(exp_i->tids.contains(al.RefID)) {
                    d=2;
                    exp_i->mapped++;
                }

                //////// HEAT MAP and ATDPes
                if(b>=0 && b<avd_whole_region) {
                    if(exp_i->regions[i]->strand) {
                        exp_i->avd_total[b]+=d;
                        if(b/avd_heat_window < exp_i->avd_matrix[i].second.size())
                            exp_i->avd_matrix[i].second[b/avd_heat_window]+=d;
                    } else {
                        exp_i->avd_total[avd_whole_region-b-1]+=d;
                        if((avd_whole_region-b-1)/avd_heat_window < exp_i->avd_matrix[i].second.size())
                            exp_i->avd_matrix[i].second[(avd_whole_region-b-1)/avd_heat_window]+=d;
                    }
                }

                ///////GENE BODY
                int _s=b1-(exp_i->regions[i]->txStart-avd_window);//
                if(_s<0) continue;
                int _idx=0;
                int gene_len=(exp_i->regions[i]->txEnd-exp_i->regions[i]->txStart);
                if(gene_len<1)  continue;
                double val=d;
                double bval=d;

                if(_s >= 0 && _s < avd_window) { //1
                    _idx=_s/r_w_b;
                    val/=r_w_b;
                } else                          //2

                    if(_s >= avd_window && _s < gene_len+avd_window ) {
                        _s -= avd_window;
                        double rat=(double)gene_len/(double)avd_bodysize;
                        _idx = _s/rat+avd_bodysize;
                        val /= rat;
                    } else                           //3

                        if(_s >=gene_len+avd_window && _s < gene_len+avd_window*2 ) {
                            _s -= (gene_len+avd_window);
                            _idx = _s/r_w_b+avd_bodysize*2;
                            val /= r_w_b;
                        } else {
                            continue;
                        }
                if(exp_i->regions[i]->strand) {
                    exp_i->body_matrix[i].second[_idx/10]=exp_i->body_matrix[i].second.at(_idx/10).toDouble()+bval;
                    exp_i->avd_body[_idx]+=val;
                } else {
                    exp_i->body_matrix[i].second[300-_idx/10-1]=exp_i->body_matrix[i].second.at(300-_idx/10-1).toDouble()+bval;
                    exp_i->avd_body[avd_bodysize*3-_idx-1]+=val;
                }

            }//trough bam reads while

        }
    }//go trough region

}

/*************
 *
 *
 *   Static functions
 *
 *
 **************/

void ATDPBasics::prn_debug(QString str,BamAlignment &al) {
    qDebug()
            <<"\n"<<str<<       al.Name.c_str()
           <<"\n IsDuplicate:"<<   al.IsDuplicate()
          <<"\n IsFailedQC"<<     al.IsFailedQC()
         <<"\n IsMaped:"<<       al.IsMapped()
        <<"\n isFirstMate:"<<   al.IsFirstMate()
       <<"\n isSecondMate:"<<  al.IsSecondMate()
      <<"\n IsMateMapped:"<<  al.IsMateMapped()
     <<"\n IsMateReverseStrand:"<< al.IsMateReverseStrand() // returns true if alignment's mate mapped to reverse strand
    <<"\n IsReverseStrand:"<<   al.IsReverseStrand()
    <<"\n IsPaired:"<<      al.IsPaired()
    <<"\n IsPrimaryAlignment:"<<al.IsPrimaryAlignment()
    <<"\n IsProperPair:"<<  al.IsProperPair()
    <<"\n AligmentFlag:"<<  al.AlignmentFlag
    <<"\n MateRefId:"<<     al.MateRefID
    <<"\n Len:"<<       al.Length
    <<"\n InsertSize:"<<    al.InsertSize
    <<"\n MatePosition:["<<al.MatePosition<<"] "
    <<"\n Position:["<<al.RefID<<":"<<al.Position+1<<"-"<<al.GetEndPosition()<<"] ";
}

