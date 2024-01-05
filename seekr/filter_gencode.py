###################################################################################################
### Description:
# filter gencode fasta file by 1: length; 2: transcript feature type and Ensemble_canonical tag 3: transcript feature type and isoform number

### Details:
# please use gencode downloaded fasta and gtf files as input, other formats are not tested
# please make sure the fasta and gtf files are from the same release and same species
# if canonical filtering is turned on, only sequences that has the feature type 'transcript' and has a tag 'Ensembl_canonical' will be kept
# if isoform filtering is turned on, only sequences that has the feature type 'transcript' and has a transcript_name ('Gm37180-201') that contains the isoform number ('201') will be kept
# isoform number support regular expression, for example, '[0-9]01' will match '101', '201', '301' up to '901'.
# if rm_dup is set to True, sequences that are exactly the same will be removed and only the first occurence will be kept
# if more than 50 transcript_ids in gtf file (after filtering by either Ensemble_canonical tag or isoform) cannot be matched to the input fasta headers, a warning message will be printed
# transcript_id in gtf file is located inside the 9th field
# transcript_id in the fasta headers is the first element of the header split by '|'
# transcript_name in gtf file is also located inside the 9th field
# these are standard formats from gencode


### Input:
# fasta_path: path to the fasta file that needs to be filtered
# gtf_path: path to the gtf file that corresponds to/matches the fasta file
# if canonical filtering is not needed, then set gtf_path=None, meaning you don't need to provide a gtf file
# len_threshold: the length threshold for filtering, if no length filtering is needed, set len_threshold=0 (default)
# canonical: whether to filter by feature type as 'transcript' and tag as 'Ensemble_canonical', if no canonical filtering is needed, set canonical=False (default)
# if canonical filtering is needed, please make sure to provide a gtf file that matches the fasta file
# if canonical filtering is turned on, only sequences that has the feature type 'transcript' and has a tag 'Ensembl_canonical' will be kept
# isoform: whether to filter by isoform number, if no isoform filtering is needed, set isoform='0' (default)
# for older version of gencode gtfs, the tag 'Ensemble_canonical' does not exist. In this case, please use isoform filtering
# accodring to gencode format, isoform must be a 3 digit number string, for example, '201', '202', '203', etc.
# isoform number can accomodate wildcard * for digit 0-9, for example, '*01' will match '101', '201', '301', etc.
# if isoform filtering is turned on, only sequences that has the feature type 'transcript' 
# and has a transcript_name ('Gm37180-201') that contains the isoform number ('201') will be kept
# please refer to Gencode data format for further details: https://www.gencodegenes.org/pages/data_format.html
# rm_dup: whether to remove duplicated sequences, rm_dup=False is the default. 
# If rm_dup=True, sequences that are exactly the same will be removed and only the first occurence will be kept
# outputname: the name of the output fasta file, default is 'test', a trailing part '.fa' will be added to the outputname to ensure fasta format


#### Output:
# fasta file (.fa) with the name of outputname.fa that is filtered by length and/or by transcript feature type and Ensemble_canonical tag and/or isoform number

### Example:
# headers, seqs = filter_gencode(fasta_path='gencode.vM33.lncRNA_transcripts.fa', 
#                                gtf_path='gencode.vM33.long_noncoding_RNAs.gtf',
#                                len_threshold=500, canonical=True, isoform='201',
#                                rm_dup=True, outputname='test')


#######################################################################################################################

from seekr.fasta_reader import Reader as seekrReader
import re

# Function to parse the 9th field of gtf 
# return the 'transcript_id' if there exists a 'tag' that has 'Ensembl_canonical' value
def get_transcript_id_with_ensembl_canonical(field):
    # Initialize transcript_id and flag for Ensembl_canonical
    transcript_id = None
    ensembl_canonical_present = False

    # Splitting the field into key-value pairs
    key_value_pairs = [kv.strip() for kv in field.split(';') if kv]
    
    for kv in key_value_pairs:
        # Splitting each pair into key and value
        key, value = kv.split(None, 1)
        value = value.strip(' "')

        # Check for transcript_id
        if key == 'transcript_id':
            transcript_id = value
        
        # Check for Ensembl_canonical tag
        if key == 'tag' and 'Ensembl_canonical' in value:
            ensembl_canonical_present = True

    # Return transcript_id if Ensembl_canonical is present, else return empty string
    return transcript_id if ensembl_canonical_present else ''

# Function to parse the 9th field of gtf 
# return the 'transcript_id' if transcript_name contains the isoform number


def get_transcript_id_with_isoform(field, isoform):
    # Initialize transcript_id and flag for Ensembl_canonical
    transcript_id = None
    isoform_match = False

    # Splitting the field into key-value pairs
    key_value_pairs = [kv.strip() for kv in field.split(';') if kv]
    
    for kv in key_value_pairs:
        # Splitting each pair into key and value
        key, value = kv.split(None, 1)
        value = value.strip(' "')

        # Check for  transcript_id
        if key == 'transcript_id':
            transcript_id = value
        
        # Check for isoform number
        if key == 'transcript_name':
            iso = value.split('-')[-1]
            # check if this iso is a 3 digit number string
            if iso.isdigit() and len(iso) == 3:
                isoform_match = bool(re.match(f'^{isoform}$', iso))

    # Return transcript_id if Ensembl_canonical is present, else return empty string
    return transcript_id if isoform_match else ''



def filter_gencode(fasta_path, gtf_path=None, len_threshold=0, canonical=False, isoform='0', rm_dup=False, outputname='test'):

    # read in fasta file
    seqs = seekrReader(fasta_path).get_seqs()
    headers=seekrReader(fasta_path).get_headers()
    # remove the '>' in the header
    headers = [i[1:] for i in headers]

    # split the headers by | and return the first element as the transcript_id and second to last element as len
    headers_list = [i.split('|') for i in headers]
    headers_tids = [i[0] for i in headers_list]
    headers_len = [int(i[-2]) for i in headers_list]
 
    # filter by canonical and/or isoform
    if canonical == True or isoform != '0':

        if gtf_path == None:
            print('Please provide a gtf file path for filtering by Ensemble_canonical tag and/or isoform number')
            return

        # read in gtf file
        with open(gtf_path, 'r') as f:
            gtfs = f.readlines()

        # filter gtf file
        # remove lines start with # which means a comment
        gtfs = [line for line in gtfs if line[0] != '#']
        # remove any leading or trailing whitespace (including newline characters). 
        # Then splits the line into a list of strings using the tab character (\t) as the delimiter
        gtfs = [line.strip().split('\t') for line in gtfs]
        # keeps only those lines where the feature type is 'transcript'
        gtfs = [line for line in gtfs if line[2] == 'transcript']

        if canonical == True: 
            # return the transcirp_id of the lines that has a tag 'Ensembl_canonical'
            tids = [get_transcript_id_with_ensembl_canonical(line[8]) for line in gtfs]
            # remove all '' from tids
            tids = [x for x in tids if x != '']

            # get the boolean of headers_tids, whether they are in tids
            # and use the positions to filter headers_tids and seqs
            # positions = [headers_tids.index(tid) for tid in tids]
            tids_set = set(tids)
            presence = [tid in tids_set for tid in headers_tids]

            # check if there are more than 50 elements in tids_set that is not in headers_tids
            # print a message and abort
            if len(tids_set - set(headers_tids)) > 50:
                print('After Ensemble_canonical tag filtering on gtf, there are more than 50 transcript_ids in gtf file that cannot be matched to the inpt fasta headers.')
                print('Please make sure the provided gtf file and fasta file are from the same release and same species.')
                print('Please use gtf and fasta files directly from gencode, other formats are not tested.')

            # Filtering headers and seq based on the presence list
            headers = [header for header, present in zip(headers, presence) if present]
            seqs = [seq for seq, present in zip(seqs, presence) if present]
            headers_len = [header_len for header_len, present in zip(headers_len, presence) if present]
            headers_tids = [tid for tid, present in zip(headers_tids, presence) if present]
            gtfs = [gtf for gtf, present in zip(gtfs, presence) if present]

        if isoform != '0':
            # return the transcirp_id of the lines that has the specified isoform number in transcript_name
            itids = [get_transcript_id_with_isoform(line[8], isoform) for line in gtfs]
            # remove all '' from itids
            itids = [x for x in itids if x != '']

            # get the boolean of headers_tids, whether they are in itids
            # and use the positions to filter headers_tids and seqs
            # positions = [headers_tids.index(tid) for tid in tids]
            itids_set = set(itids)
            ipresence = [tid in itids_set for tid in headers_tids]

            # check if there are more than 50 elements in itids_set that is not in headers_tids
            # print a message and abort
            if len(itids_set - set(headers_tids)) > 50:
                print('After isoform filtering on gtf, there are more than 50 transcript_ids in gtf file that cannot be matched to the input fasta headers.')
                print('Please make sure the provided gtf file and fasta file are from the same release and same species.')
                print('Please use gtf and fasta files directly from gencode, other formats are not tested.')
            
            # Filtering headers and seq based on the ipresence list
            headers = [header for header, ipresent in zip(headers, ipresence) if ipresent]
            seqs = [seq for seq, ipresent in zip(seqs, ipresence) if ipresent]
            headers_len = [header_len for header_len, ipresent in zip(headers_len, ipresence) if ipresent]



    # filter by length
    if len_threshold > 0:
        seqs = [seq for seq, header_len in zip(seqs, headers_len) if header_len >= len_threshold]
        headers = [header for header, header_len in zip(headers, headers_len) if header_len >= len_threshold]

    # remove duplicated sequences and keep the corresponding headers
    if rm_dup == True:
        seen_seqs = set()
        seqs_uni = []
        headers_uni = []

        for seq, header in zip(seqs, headers):
            if seq not in seen_seqs:
                seen_seqs.add(seq)
                seqs_uni.append(seq)
                headers_uni.append(header)

        seqs = seqs_uni
        headers = headers_uni
        
   
    # write out filtered fasta file
    with open(f'{outputname}.fa', 'w') as f:
        for i in range(len(headers)):
            f.write('>'+headers[i]+'\n'+seqs[i]+'\n')
            
    return headers, seqs

