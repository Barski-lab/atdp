/****************************************************************************
**
** Copyright (C) 2011-2014 Andrey Kartashov .
** All rights reserved.
** Contact: Andrey Kartashov (porter@porter.st)
**
** This file is part of the averagedensity module of the genome-tools.
**
** GNU Lesser General Public License Usage
** This file may be used under the terms of the GNU Lesser General Public
** License version 2.1 as published by the Free Software Foundation and
** appearing in the file LICENSE.LGPL included in the packaging of this
** file. Please review the following information to ensure the GNU Lesser
** General Public License version 2.1 requirements will be met:
** http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html.
**
** Other Usage
** Alternatively, this file may be used in accordance with the terms and
** conditions contained in a signed written agreement between you and Andrey Kartashov.
**
****************************************************************************/

#include "atdp.hpp"

class ATDPThreads: public QRunnable
{

    public:

        EXPERIMENT_INFO* exp_i;
        ATDPBasics* atdpb;

        ATDPThreads(EXPERIMENT_INFO* e):
            exp_i(e) {
            atdpb = new ATDPBasics(exp_i);
            this->setAutoDelete(true);
        };

        void run(void){
            atdpb->RegionsProcessing();
        }

};

//-------------------------------------------------------------
//-------------------------------------------------------------
ATDP::ATDP(QObject* parent):
    QObject(parent)
{        
    avd_window=gArgs().getArgs("avd_window").toInt();
    avd_whole_region=avd_window*2+1;
    avd_heat_window=gArgs().getArgs("avd_heat_window").toInt();
    twicechr=gArgs().getArgs("sam_twicechr").toString();
    ignorechr=gArgs().getArgs("sam_ignorechr").toString();
    avd_bodysize=1000;
}

void ATDP::start() {
    qDebug()<<"start";

    //Prepare experiments to proccess
    this->getRecordInfo();

    QThreadPool *t_pool=QThreadPool::globalInstance();

    //foreach trough experiments
    foreach(QString key,experiment_info.keys()){
        t_pool->start(new ATDPThreads(&experiment_info[key]));
    }

    if(t_pool->activeThreadCount()!=0) {
        qDebug()<<"waiting threads";
        t_pool->waitForDone();
        qDebug()<<"threads done";
    }

    if ( !export_to_file( gArgs().getArgs("out").toString().toStdString() ) ){
        throw "Error export results";
    }


//    } else { // Advanced Analyses

//        QJsonArray columns_name;
//        for(int col=-avd_window;col<avd_window;col+=avd_heat_window)
//            columns_name.append(QString(""));
//        columns_name[0]=QString("-%1k").arg(avd_window/1000);
//        columns_name[columns_name.count()-1]=QString("+%1k").arg(avd_window/1000);
//        columns_name[columns_name.count()/2]=QString("TSS");

//        //QJsonArray data_array;
//        QList<QJsonObject*> data_array_obj;

//        foreach(QString key,experiment_info.keys()){
//            exp_i=&experiment_info[key];
//            QJsonObject *s_data = new QJsonObject();
//            QJsonObject &data = (*s_data);
//            data["cols"]=columns_name;
//            data["tbl1_id"]=exp_i->tbl1_id;
//            data["tbl2_id"]=exp_i->tbl2_id;
//            data["tbl1_name"]=exp_i->tbl1_name;
//            data["tbl2_name"]=exp_i->tbl2_name;
//            data["pltname"]=exp_i->plotname;
//            data["mapped"]=exp_i->mapped;

//            QJsonArray matrix;
//            QJsonArray body_matrix;
//            QJsonArray rpkm_matrix;
//            QJsonArray rows;
//            QJsonArray glength;
//            QJsonArray gene_body;
//            /*
//             * Avd gene body
//             */
//            QVector<double> storage;
//            for(int w=0; w < exp_i->avd_body.size(); w++)
//                storage<<(exp_i->avd_body.at(w)/exp_i->mapped)/exp_i->regions.size();
//            storage=Math::smooth<double>(storage,gArgs().getArgs("avd_bsmooth").toInt());
//            for(int w=0; w<storage.size(); w++)
//                gene_body.append(storage.at(w));
//            /*
//             * Gene body matrix
//             */
//            for(int j=0; j<exp_i->body_matrix.size();j++) {
//                body_matrix.append(exp_i->body_matrix[j].second);
//            }



//            /*
//             *  AVD HEAT
//             */
//            QList<quint64> max;
//            for(int j=0; j<exp_i->avd_matrix.size();j++) {
//                QJsonArray row;
//                quint64 max_line=0;
//                for(int c=0; c<exp_i->avd_matrix[j].second.size();c++) {
//                    max_line=qMax<double>(exp_i->avd_matrix[j].second[c],max_line);//maximum value in a row
//                    row.append(exp_i->avd_matrix[j].second[c]);
//                }

//                max<<max_line;
//                matrix.append(row);
//                rows.append(exp_i->rpkm_matrix[j].first->gene_id);
//                glength.append(exp_i->rpkm_matrix[j].first->txEnd-exp_i->rpkm_matrix[j].first->txStart);
//                rpkm_matrix.append(exp_i->rpkm_matrix[j].second);
//            }

//            qSort(max);
//            float quan=qCeil(max.size()*0.9);
//            data["rows"]=rows;
//            data["glengths"]=glength;
//            double maxx=(double)max.at(quan);
//            if(quan>1) {
//                maxx+=max.at(quan-1);
//                maxx/=2;
//            }
//            data["max"]=maxx;
//            data["array"]=matrix;
//            data["bodyarray"]=body_matrix;
//            data["rpkmarray"]=rpkm_matrix;
//            data["rpkmcols"]=exp_i->rpkmnames;
//            data["genebody"]=gene_body;
//            data_array_obj<<s_data;
//        }//foreach trough experiments

//        QFile outFile;
//        QString out="{\"message\":\"Data populated\",\"success\":true,\"total\":";
//        outFile.open(stdout,QIODevice::WriteOnly);
//        out+=QString("%1,\"data\":[").arg(experiment_info.size());
//        outFile.write(out.toLocal8Bit());

//        for(int ii=0;ii<data_array_obj.size(); ii++) {
//            outFile.write(QJsonDocument(*(data_array_obj.at(ii))).toJson(QJsonDocument::Compact));
//            if(ii<data_array_obj.size()-1)
//                outFile.write(",\n");
//        }
//        outFile.write("]}\n");
//        outFile.close();
//    }


    qDebug()<<"end";

    emit finished();
}


void ATDP::getRecordInfo() {
    if(gArgs().getArgs("in").toString().isEmpty()) {
        throw "Set input BAM file";
    }

    EXPERIMENT_INFO *ei = new EXPERIMENT_INFO();
    ei->fragmentsize=gArgs().getArgs("fragmentsize").toInt();
    ei->filepath=gArgs().getArgs("in").toString();
    ei->avd_total.resize(avd_whole_region);
    ei->avd_total.fill(0,avd_whole_region);
    ei->avd_body.resize(avd_bodysize*3);
    ei->avd_body.fill(0,avd_bodysize*3);
    ei->mapped=gArgs().getArgs("mappedreads").toInt();
    experiment_info.insert(gArgs().getArgs("in").toString(),*ei);
}

/*
 * designed for atdp advanced analysis takes list pairs <tbl1_id, tbl2_id> (<chip,rna>)
 */


//void ATDP::getRecordsInfo() {
//    QSqlQuery q;

//    q.prepare("select a.pltname,tbl1_id,tbl2_id,l.fragmentsize,l.tagsmapped,g2.tableName,g.tableName,g2.name,g.name "
//              " from atdp a,genelist g, genelist g2, labdata l where a.genelist_id=? "
//              " and g.labdata_id=l.id and a.tbl1_id=g.id and a.tbl2_id=g2.id");
//    q.bindValue(0, gArgs().getArgs("avd_guid").toString());
//    if(!q.exec()) {
//        qDebug()<<"Query error info: "<<q.lastError().text();
//        throw "Error query to DB";
//    }
//    while(q.next()) {
//        EXPERIMENT_INFO *ei = new EXPERIMENT_INFO();
//        ei->source=q.value(5).toString();
//        ei->db=gSettings().getValue("experimentsdb");
//        ei->fragmentsize=q.value(3).toInt();
//        ei->mapped=q.value(4).toInt();
//        ei->filepath=q.value(6).toString()+"/"+q.value(6).toString()+".bam";
//        ei->avd_total.resize(avd_whole_region);
//        ei->avd_total.fill(0,avd_whole_region);
//        ei->avd_body.resize(avd_bodysize*3);
//        ei->avd_body.fill(0,avd_bodysize*3);
//        ei->plotname=q.value(0).toString();
//        ei->tbl1_id=q.value(1).toString();
//        ei->tbl2_id=q.value(2).toString();
//        ei->tbl1_name=q.value(8).toString();
//        ei->tbl2_name=q.value(7).toString();
//        experiment_info.insert(q.value(2).toString()+"="+q.value(1).toString(),*ei);
//    }
//}

void ATDP::print(ostream& output_stream){
    EXPERIMENT_INFO* exp_i;
    exp_i=&experiment_info[experiment_info.keys().at(0)];
    QVector<double> storage;
    for(int w=0; w< exp_i->avd_total.size(); w++)
        storage<<(exp_i->avd_total.at(w)/exp_i->mapped)/exp_i->regions.size();
    storage=Math::smooth<double>(storage,gArgs().getArgs("avd_smooth").toInt());
    int rows=storage.size();
    output_stream // header line
            << "X" << "\t"
            << "Y" << endl;
    for(int i=0; i<rows;i++) {
        output_stream << (int)(i-rows/2) << "\t";
        output_stream << std::setprecision (15) << storage.at(i) << endl;
    }
}

bool ATDP::export_to_file(const string & output_filename){
    ofstream output_stream (output_filename.c_str());
    if (output_stream.is_open())
    {
        print(output_stream);
        output_stream.close();
        return true;
    }
    return false;
}



ATDP::~ATDP()
{
}


