#include "bam_reader_util.hpp"




// Check if current bam file is indexed (and that index data is loaded into program)
bool make_index (BamReader & bam_reader){
    if (not bam_reader.HasIndex()){
        cerr << "Current BAM file isn't indexed" << endl;
        cerr << "Trying to find index files in the same directory" << endl;
        // Trying to load index data from the filesystem
        // If BamIndex::STANDARD we are looking for BAI file
        if (bam_reader.LocateIndex(BamIndex::STANDARD)){
            cerr << "Located and loaded index file from disk" << endl;
        } else {
            cerr << "Couldn't locate index file" << endl;
            cerr << "Trying to create the new one" << endl;
            // Trying to create index data ourself
            // BamIndex::STANDARD - we are trying to create BAI file
            if (not bam_reader.CreateIndex(BamIndex::STANDARD)){
                cerr << "Cannot create index for current bam file. Exiting" << endl;
                return false;
            };
            cerr << "Index file for current BAM file is succesfully created" << endl;
        }
    }
    return true;
}

void get_bam_info(BamReader & bam_reader, BamGeneralInfo & bam_general_info){

    BamAlignment al;
    while (bam_reader.GetNextAlignment(al)){
        flag_check (al, bam_general_info);
    }
    bam_general_info.aligned = bam_general_info.total - bam_general_info.not_aligned;
    bam_reader.Rewind();

    cout << "BAM alignment statistics:" << endl;
    cout << "   Total: " << bam_general_info.total << endl;
    cout << "   Aligned: " << bam_general_info.aligned << endl;
    cout << "   Not aligned: " << bam_general_info.not_aligned << endl;
}


bool flag_check (BamAlignment & al, BamGeneralInfo & bam_general_info){
    bam_general_info.total++;
    if( al.IsMapped() && al.IsPrimaryAlignment() && (!al.IsDuplicate()) ) {
        if (al.IsPaired()){
            if(al.IsProperPair() && al.IsMateMapped() ) {
                // The only possible situation when it may happen is when read was mapped on the wrong strand
                if( (!(al.IsReverseStrand() ^ al.IsMateReverseStrand())           ) ||  // both mates came from the same strand. it shouldn't be like this
                    ( (al.Position < al.MatePosition) && al.IsReverseStrand()     ) ||
                    ( (al.MatePosition < al.Position) && al.IsMateReverseStrand() )
                        ) {
                    bam_general_info.not_aligned++;
                    return false;
                }
            } else {
                bam_general_info.not_aligned++;
                return false;
            }
        }
    } else {
        bam_general_info.not_aligned++;
        return false;
    }
    return true;
}
