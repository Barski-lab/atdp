#ifndef BAM_READER_UTIL_HPP
#define BAM_READER_UTIL_HPP

#endif // BAM_READER_UTIL_HPP

#include <api/BamMultiReader.h>
#include <api/BamWriter.h>

using namespace BamTools;
using namespace std;


struct BamGeneralInfo {
    long total;
    long not_aligned;
    long aligned;
    BamGeneralInfo ():
            total (0),
            not_aligned (0),
            aligned (0)
    {
    }
    void operator = (const BamGeneralInfo &other ) {
        total = other.total;
        not_aligned = other.not_aligned;
        aligned = other.aligned;
    }
};

bool make_index (BamReader & bam_reader);
void get_bam_info(BamReader & bam_reader, BamGeneralInfo & bam_general_info);
bool flag_check (BamAlignment & al, BamGeneralInfo & bam_general_info);
