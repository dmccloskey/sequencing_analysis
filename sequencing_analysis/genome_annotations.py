from io_utilities.base_importData import base_importData
from io_utilities.base_exportData import base_exportData
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqFeature
from Bio.SeqFeature import FeatureLocation

class genome_annotations():
    '''class for genome annotations'''
    def __init__(self,annotation_I=None,annotation_ref_I=None,geneReference_I=None,sequence_I=None,sequence_ref_I=None,
                 IS_sequences_I=None,IS_sequences_ref_I=None):
        
        if annotation_I and annotation_ref_I:
            self.annotation = self.set_annotation(annotation_I,annotation_ref_I);
        else:
            self.annotation = None;
        if sequence_I and sequence_ref_I:
            self.sequence = self.set_sequence(sequence_I,sequence_ref_I);
        else:
            self.sequence = None;
        if IS_sequences_I and IS_sequences_ref_I:
            self.IS_sequences = self.set_sequences(IS_sequences_I,IS_sequences_ref_I);
        else:
            self.IS_sequences = {};
        '''
        geneReference
        d = {};
        d['biologicalmaterial_id'],
        d['ordered_locus_name'],
        d['ordered_locus_name2'],
        d['swissprot_entry_name'],
        d['ac'],
        d['ecogene_accession_number'],
        d['gene_name'])'''
        if geneReference_I:
            self.geneReference = self.import_biologicalMaterialGeneReferences(geneReference_I);
        else:
            self.geneReference = [];

    def import_biologicalMaterialGeneReferences(self, filename):
        '''import biological material gene references for annotations'''
        io = base_importData();
        io.read_csv(filename);
        return io.data;

    def set_annotation(self,ref_genome_I,ref_I='genbank'):
        """ref_genome_I = reference genome to use for the annotation
        ref_I = reference database, default, 'genbank'
        """
        return SeqIO.read(ref_genome_I,ref_I);

    def set_sequence(self,ref_genome_I,ref_I='fasta'):
        """ref_genome_I = reference genome sequence
        ref_I = reference database, default, 'fasta'
        """
        return SeqIO.read(ref_genome_I,ref_I);

    def set_sequences(self,sequences_I,ref_I='fasta'):
        """sequences_I = list of sequences
        ref_I = reference database, default, 'fasta'
        """
        '''TODO:
        from Bio import SeqIO
            for index, record in enumerate(SeqIO.parse(open("ls_orchid.gbk"), "genbank")):
                print "index %i, ID = %s, length %i, with %i features" \
                      % (index, record.id, len(record.seq), len(record.features))
        '''
        records = {};
        for record in SeqIO.parse(open(sequences_I),ref_I):
            records[record.id] = record;
        return records;
        
    def _find_genesFromMutationPosition(self,mutation_position_I,annotation_I=None):
        '''find genes at the position or closest to the position given the reference genome
        input:
        mutation_position_I = mutation position [int]
        annotation = genbank record for the annotation [SeqRecord]
        '''
        if not annotation_I:
            annotation_I = self.annotation;
        snp_records = {};
        snp_records['gene'] = []
        snp_records['db_xref'] = []
        snp_records['locus_tag'] = []
        snp_records['EC_number'] = []
        snp_records['product'] = []
        snp_records['location'] = []
        # find features in the coding region of the genome that bracket the mutation position
        for feature_cnt,feature in enumerate(annotation_I.features):
            if mutation_position_I in feature and feature.type == 'gene':
                snp_records['gene'] = feature.qualifiers.get('gene')
                snp_records['db_xref'] = feature.qualifiers.get('db_xref')
                snp_records['locus_tag'] = feature.qualifiers.get('locus_tag')
            elif mutation_position_I in feature and feature.type == 'CDS':
                if feature.qualifiers.get('EC_number'):snp_records['EC_number'] = feature.qualifiers.get('EC_number')
                else:snp_records['EC_number'] = [None];
                if feature.qualifiers.get('product'):snp_records['product'] = feature.qualifiers.get('product')
                else:snp_records['product'] = [None];
                snp_records['location'] = ['coding'];
            elif mutation_position_I in feature and feature.type == 'repeat_region':
                snp_records['location'] = feature.qualifiers.get('note')
            elif mutation_position_I in feature and feature.type == 'mobile_element':
                snp_records['location'] = feature.qualifiers.get('mobile_element_type')
            elif mutation_position_I in feature and feature.type == 'misc_feature':
                snp_records['location'] = feature.qualifiers.get('note')
            elif mutation_position_I in feature and feature.type == 'mat_peptide':
                snp_records['gene'] = feature.qualifiers.get('gene')
                snp_records['locus_tag'] = feature.qualifiers.get('locus_tag')
                #snp_records['location'] = feature.qualifiers.get('note')
                if feature.qualifiers.get('EC_number'):snp_records['EC_number'] = feature.qualifiers.get('EC_number')
                else:snp_records['EC_number'] = [None];
                if feature.qualifiers.get('product'):snp_records['product'] = feature.qualifiers.get('product')
                else:snp_records['product'] = [None];
                snp_records['location'] = ['coding'];
            elif mutation_position_I in feature and feature.type == 'tRNA':
                snp_records['gene'] = feature.qualifiers.get('gene')
                snp_records['locus_tag'] = feature.qualifiers.get('locus_tag')
                #snp_records['location'] = feature.qualifiers.get('note')
                if feature.qualifiers.get('EC_number'):snp_records['EC_number'] = feature.qualifiers.get('EC_number')
                else:snp_records['EC_number'] = [None];
                if feature.qualifiers.get('product'):snp_records['product'] = feature.qualifiers.get('product')
                else:snp_records['product'] = [None];
                snp_records['location'] = ['coding'];
            elif mutation_position_I in feature and feature.type == 'rRNA':
                snp_records['gene'] = feature.qualifiers.get('gene')
                snp_records['db_xref'] = feature.qualifiers.get('db_xref')
                snp_records['locus_tag'] = feature.qualifiers.get('locus_tag')
                #snp_records['location'] = feature.qualifiers.get('note')
                if feature.qualifiers.get('EC_number'):snp_records['EC_number'] = feature.qualifiers.get('EC_number')
                else:snp_records['EC_number'] = [None];
                if feature.qualifiers.get('product'):snp_records['product'] = feature.qualifiers.get('product')
                else:snp_records['product'] = [None];
                snp_records['location'] = ['coding'];
            elif mutation_position_I in feature and feature.type == 'ncRNA':
                snp_records['gene'] = feature.qualifiers.get('gene')
                snp_records['db_xref'] = feature.qualifiers.get('db_xref')
                snp_records['locus_tag'] = feature.qualifiers.get('locus_tag')
                #snp_records['location'] = feature.qualifiers.get('note')
                if feature.qualifiers.get('EC_number'):snp_records['EC_number'] = feature.qualifiers.get('EC_number')
                else:snp_records['EC_number'] = [None];
                if feature.qualifiers.get('product'):snp_records['product'] = feature.qualifiers.get('product')
                else:snp_records['product'] = [None];
                snp_records['location'] = ['coding'];
            elif mutation_position_I in feature and feature.type != 'source':
                print(feature)
        if not snp_records['location']:
            # no features in the coding region were found that bracket the mutation
            # find features before and after the mutation position
            start_prev = 0;
            stop_prev = 0;
            inter1_start = None;
            inter1_stop = None;
            inter2_start = None;
            inter2_stop = None;
            # pass 1: locate the start and stop positions of the features before and after the mutation
            for feature_cnt,feature in enumerate(annotation_I.features):
                start = feature.location.start.position
                stop = feature.location.end.position
                if mutation_position_I > stop_prev and mutation_position_I < start:
                    inter1_start = start_prev;
                    inter1_stop = stop_prev;
                    inter2_start = start;
                    inter2_stop = stop
                    break;
                start_prev = start;
                stop_prev = stop;
            if not inter1_start:
                # the end of the genome was reached without finding features both before and after the mutation
                # record the last entry
                inter1_start = start_prev;
                inter1_stop = stop_prev;
                inter2_start = start;
                inter2_stop = stop
            # pass 2: record features before and after the mutation
            for feature_cnt,feature in enumerate(annotation_I.features):
                start = feature.location.start.position
                stop = feature.location.end.position
                if (inter1_start == start and inter1_stop == stop) or (inter2_start == start and inter2_stop == stop):
                    if feature.type == 'gene':
                        snp_records['gene'] += feature.qualifiers.get('gene')
                        snp_records['db_xref'] += feature.qualifiers.get('db_xref')
                        snp_records['locus_tag'] += feature.qualifiers.get('locus_tag')
                    if feature.type == 'CDS':
                        if feature.qualifiers.get('EC_number'):snp_records['EC_number'] += feature.qualifiers.get('EC_number')
                        else:snp_records['EC_number'] += [None]
                        if feature.qualifiers.get('product'):snp_records['product'] += feature.qualifiers.get('product')
                        else:snp_records['product'] += [None]
            for gene in snp_records['gene']:
                snp_records['location'] += ['intergenic']
        return snp_records;
    def _generate_httplink2gene_ecogene(self,ecogene_I):
        '''Generate link to ecocyc using the ecogene accession number'''
        ecogene_httplink = 'http://ecocyc.org/ECOLI/NEW-IMAGE?type=GENE&object='+ecogene_I;
        return ecogene_httplink
    def _find_genesInRegion(self,start_I,stop_I,annotation_I=None):
        '''find genes in the start and stop region of the genome
        INPUT:
        start_I = mutation start position [int]
        stop_I = mutation stop position [int]
        annotation_I = genbank record annotations [SeqRecord]
        '''
        if not annotation_I:
            annotation_I = self.annotation;
        data_O = [];
        #extract all features within the start and stop region
        features = [f for f in annotation_I.features if start_I <= f.location.start.position and stop_I <= f.location.end.position]
        for feature_cnt,feature in enumerate(features):
            # NOTE:
            # there are two records for each gene: one of type "gene" and another of type "CDS"
            # sometimes there is also a third of type "mat_peptide" which has a different start/stop position
            # this algorithm will combine the records for types "gene" and "CDS" as a single row
            # and add mat_peptide as a second row (if desired)
            # initialize variables
            if feature_cnt == 0:
                feature_start_pos = feature.location.start.position;
                feature_stop_pos = feature.location.end.position;
                snp_records = {};
                snp_records['gene'] = []
                snp_records['db_xref'] = []
                snp_records['locus_tag'] = []
                snp_records['EC_number'] = []
                snp_records['product'] = []
                snp_records['location'] = []
                snp_records['start'] = None
                snp_records['stop'] = None
                snp_records['type'] = []
            # add complete snp_record to data_O
            if feature_start_pos != feature.location.start.position or feature_stop_pos != feature.location.end.position:
                #snp_records['start'] should be in every record
                data_O.append(snp_records);
                feature_start_pos = feature.location.start.position;
                feature_stop_pos = feature.location.end.position;
                snp_records = {};
                snp_records['gene'] = []
                snp_records['db_xref'] = []
                snp_records['locus_tag'] = []
                snp_records['EC_number'] = []
                snp_records['product'] = []
                snp_records['location'] = []
                snp_records['start'] = None
                snp_records['stop'] = None
                snp_records['type'] = []
            # fill in snp_record (may require multiple passes through features)
            if feature.type == 'gene':
                snp_records['gene'] = feature.qualifiers.get('gene')
                snp_records['db_xref'] = feature.qualifiers.get('db_xref')
                snp_records['locus_tag'] = feature.qualifiers.get('locus_tag')
                snp_records['start'] = feature.location.start.position;
                snp_records['stop'] = feature.location.end.position;
                snp_records['type'].append(feature.type)
            elif feature.type == 'CDS': 
                if feature.qualifiers.get('EC_number'):snp_records['EC_number'] = feature.qualifiers.get('EC_number')
                else:snp_records['EC_number'] = [None];
                if feature.qualifiers.get('product'):snp_records['product'] = feature.qualifiers.get('product')
                else:snp_records['product'] = [None];
                snp_records['location'] = ['coding'];
                snp_records['start'] = feature.location.start.position;
                snp_records['stop'] = feature.location.end.position;
                snp_records['type'].append(feature.type)
            elif feature.type == 'repeat_region':
                snp_records['location'] = feature.qualifiers.get('note')
                snp_records['start'] = feature.location.start.position;
                snp_records['stop'] = feature.location.end.position;
                snp_records['type'].append(feature.type)
            elif feature.type == 'mobile_element':
                snp_records['location'] = feature.qualifiers.get('mobile_element_type')
                snp_records['start'] = feature.location.start.position;
                snp_records['stop'] = feature.location.end.position;
                snp_records['type'].append(feature.type)
            elif feature.type == 'misc_feature':
                snp_records['location'] = feature.qualifiers.get('note')
                snp_records['start'] = feature.location.start.position;
                snp_records['stop'] = feature.location.end.position;
                snp_records['type'].append(feature.type)
            elif feature.type == 'mat_peptide':
                snp_records['gene'] = feature.qualifiers.get('gene')
                snp_records['locus_tag'] = feature.qualifiers.get('locus_tag')
                #snp_records['location'] = feature.qualifiers.get('note')
                if feature.qualifiers.get('EC_number'):snp_records['EC_number'] = feature.qualifiers.get('EC_number')
                else:snp_records['EC_number'] = [None];
                if feature.qualifiers.get('product'):snp_records['product'] = feature.qualifiers.get('product')
                else:snp_records['product'] = [None];
                snp_records['location'] = ['coding'];
                snp_records['start'] = feature.location.start.position;
                snp_records['stop'] = feature.location.end.position;
                snp_records['type'].append(feature.type)
            elif feature.type == 'tRNA':
                snp_records['gene'] = feature.qualifiers.get('gene')
                snp_records['locus_tag'] = feature.qualifiers.get('locus_tag')
                #snp_records['location'] = feature.qualifiers.get('note')
                if feature.qualifiers.get('EC_number'):snp_records['EC_number'] = feature.qualifiers.get('EC_number')
                else:snp_records['EC_number'] = [None];
                if feature.qualifiers.get('product'):snp_records['product'] = feature.qualifiers.get('product')
                else:snp_records['product'] = [None];
                snp_records['location'] = ['coding'];
                snp_records['start'] = feature.location.start.position;
                snp_records['stop'] = feature.location.end.position;
                snp_records['type'].append(feature.type)
            elif feature.type == 'rRNA':
                snp_records['gene'] = feature.qualifiers.get('gene')
                snp_records['db_xref'] = feature.qualifiers.get('db_xref')
                snp_records['locus_tag'] = feature.qualifiers.get('locus_tag')
                #snp_records['location'] = feature.qualifiers.get('note')
                if feature.qualifiers.get('EC_number'):snp_records['EC_number'] = feature.qualifiers.get('EC_number')
                else:snp_records['EC_number'] = [None];
                if feature.qualifiers.get('product'):snp_records['product'] = feature.qualifiers.get('product')
                else:snp_records['product'] = [None];
                snp_records['location'] = ['coding'];
                snp_records['start'] = feature.location.start.position;
                snp_records['stop'] = feature.location.end.position;
                snp_records['type'].append(feature.type)
            elif feature.type == 'ncRNA':
                snp_records['gene'] = feature.qualifiers.get('gene')
                snp_records['db_xref'] = feature.qualifiers.get('db_xref')
                snp_records['locus_tag'] = feature.qualifiers.get('locus_tag')
                #snp_records['location'] = feature.qualifiers.get('note')
                if feature.qualifiers.get('EC_number'):snp_records['EC_number'] = feature.qualifiers.get('EC_number')
                else:snp_records['EC_number'] = [None];
                if feature.qualifiers.get('product'):snp_records['product'] = feature.qualifiers.get('product')
                else:snp_records['product'] = [None];
                snp_records['location'] = ['coding'];
                snp_records['start'] = feature.location.start.position;
                snp_records['stop'] = feature.location.end.position;
                snp_records['type'].append(feature.type)
            elif feature.type != 'source':
                print(feature);
            # add the final snp_record to data_O
            if feature_cnt == len(features)-1:
                data_O.append(snp_records);

        return data_O;
    def _get_ecogenesByBiologicalmaterialIDAndOrderedLocusName(self,biologicalmaterial_id_I,ordered_locus_name_I):
        """return the ecogenes that match the biologicalmaterial_id and ordered_locus_name"""
        ecogenes = [];
        for d in self.geneReference:
            if d['biologicalmaterial_id'] ==  biologicalmaterial_id_I and d['ordered_locus_name'] == ordered_locus_name_I:
                ecogenes.append(d);
        return ecogenes;

    def convert_Seq2String(self,Seq_I):
        '''convert Seq object to string'''
        seq_str = str(Seq_I.lower());
        return seq_str;

    def convert_Seq2List(self,Seq_I):
        '''convert Seq object to list'''
        seq_list = list(Seq_I.lower());
        return seq_list;

    def make_SeqFromString(self,sequence_str_I):
        '''make a Seq object from a string'''

        sequence_O = Seq(sequence_str_I, generic_dna);
        return sequence_O;

    def check_stopCodon(self,peptide_I):
        '''Check for a stop codon at the end of the peptide sequence'''
        has_stop_codon = False;
        if peptide_I and peptide_I[-1]=='*':
            has_stop_codon = True;
        elif peptide_I and peptide_I[-1]!='*':
            print('stop codon not found');
        elif not peptide_I: 
            print('no peptide');
        return has_stop_codon;

    def transcribeAndTranslate_feature(self,sequence_I,feature_I,table_I = "Standard"):
        '''Transcribe and translate a feature from a reference sequence
        INPUT:
        sequence_I = sequence object
        feature_I = feature object
        table_I = translation table (e.g. Bacterial), default is Standard
        OUTPUT:
        feature_sequence = sequence object of the dna sequence
        rna = sequence object of the rna sequence
        peptide = sequence object of the peptide sequence
        stop_codon_found = Boolean, True if a stop codon was found
        '''

        feature_sequence = feature_I.extract(sequence_I);
        ## check for strand orientation not needed:
        #if feature_I.location.strand != -1:
        #    rna = feature_sequence.reverse_complement().transcribe();
        #else:
        #    rna = feature_sequence.transcribe();
        rna = feature_sequence.transcribe();
        peptide = rna.translate(table=table_I,to_stop=False);
        has_stop_codon = self.check_stopCodon(peptide);
        peptide = rna.translate(table=table_I,to_stop=True);

        return feature_sequence, rna, peptide, has_stop_codon;

    def mutate_sequence(self,sequence_I,mutation_I,feature_I):
        '''mutate a sequence
        INPUT:
        sequence_I = sequence object
        mutation_I = {} mutation information
        feature_I = feature object
        OUTPUT:
        sequence_O = sequence object
        mutation_O = {} updated mutation information
        feature_O = feature object
        '''
        sequence_O_str = list(sequence_I.lower());
        feature_O = None;
        mutation_O = mutation_I;

        #parse the mutation dict
        if 'position' in mutation_I: mutation_position_I = mutation_I['position'];
        if 'type' in mutation_I: mutation_type_I = mutation_I['type'];
        if 'new_seq' in mutation_I: new_seq_I = mutation_I['new_seq'].lower();
        if 'size' in mutation_I: size_I = mutation_I['size'];
        elif 'new_seq' in mutation_I: size_I = len(new_seq_I);
        else: size_I = 0;
        if 'new_copy_number' in mutation_I: new_copy_number_I = mutation_I['new_copy_number'];

        if mutation_type_I=='SNP':
            sequence_O_str[mutation_position_I-1]=new_seq_I;
            feature_O = feature_I;
        elif mutation_type_I=='DEL':
            sequence_O_str=sequence_O_str[:mutation_position_I-1]+sequence_O_str[mutation_position_I-1+size_I:];
            feature_O = self.make_SeqFeature(feature_I.location.start,
                                             feature_I.location.end - size_I,
                                             feature_I.strand,
                                             feature_I.type); # adjust the feature length
            mutation_O['feature_position'] = self.convert_genomePosition2FeaturePosition(mutation_position_I,feature_I);
        elif mutation_type_I=='INS':
            sequence_O_str = sequence_O_str[:mutation_position_I-1] + list(new_seq_I) + sequence_O_str[mutation_position_I-1:];
            feature_O = self.make_SeqFeature(feature_I.location.start,
                                             feature_I.location.end + len(new_seq_I),
                                             feature_I.strand,
                                             feature_I.type); # adjust the feature length
            mutation_O['feature_position'] = self.convert_genomePosition2FeaturePosition(mutation_position_I,feature_I);
        elif mutation_type_I=='SUB':
            sequence_O_str[mutation_position_I-1:mutation_position_I-1+size_I]=list(new_seq_I);
            end = feature_I.location.end + len(new_seq_I) - size_I; # deletion + insertion
            feature_O = self.make_SeqFeature(feature_I.location.start,
                                             end,
                                             feature_I.strand,
                                             feature_I.type); # adjust the feature length
            mutation_O['feature_position'] = self.convert_genomePosition2FeaturePosition(mutation_position_I,feature_I);
        elif mutation_type_I=='AMP':
            amp_seq = sequence_O_str[mutation_position_I-1:mutation_position_I-1+size_I];
            new_seq_I = [];
            for n in range(new_copy_number_I):
                new_seq_I.extend(amp_seq);
            size_I = len(new_seq_I);
            sequence_O_str = sequence_O_str[:mutation_position_I-1] + new_seq_I + sequence_O_str[mutation_position_I-1+size_I:];
            feature_O = self.make_SeqFeature(feature_I.location.start,
                                             feature_I.location.end + len(new_seq_I),
                                             feature_I.strand,
                                             feature_I.type); # adjust the feature length
            mutation_O['feature_position'] = self.convert_genomePosition2FeaturePosition(mutation_position_I,feature_I);
            mutation_O['new_feature_sequence'] = ''.join(new_seq_I);
        elif mutation_type_I=='MOB':
            is_name = mutation_I['repeat_name'];
            duplication_size = mutation_I['duplication_size'];
            is_sequence = self.get_ISsequence(is_name);
            if is_sequence:
                is_sequence_list = self.convert_Seq2List(is_sequence);
            else:
                return None,feature_O;
            new_seq_I = [];
            for n in range(duplication_size):
                new_seq_I.extend(is_sequence_list);
            size_I = len(new_seq_I);
            sequence_O_str = sequence_O_str[:mutation_position_I-1] + new_seq_I + sequence_O_str[mutation_position_I-1:];
            feature_O = self.make_SeqFeature(feature_I.location.start,
                                                feature_I.location.end + len(new_seq_I),
                                                feature_I.strand,
                                                feature_I.type); # adjust the feature length
            mutation_O['feature_position'] = self.convert_genomePosition2FeaturePosition(mutation_position_I,feature_I);
            mutation_O['new_feature_sequence'] = ''.join(new_seq_I);
        else:
            print('mutation type ' +mutation_type_I+' not yet supported!');
            return None,feature_O,mutation_O;

        sequence_O_str = ''.join(sequence_O_str);
        sequence_O = self.make_SeqFromString(sequence_O_str);
        return sequence_O,feature_O,mutation_O;

    def get_ISsequence(self,is_name_I):
        '''Get the IS sequence'''
        if is_name_I in self.IS_sequences:
            is_sequence = self.IS_sequences[is_name_I]
            return is_sequence.seq;
        else:
            print("IS sequence " + is_name_I + " not found.");
            return None;

    def make_SeqFeature(self,start_pos_I,stop_pos_I,strand_I,type_I,
                        location_operator_I='',qualifiers_I=None,sub_features_I=None,
                        ref_I = None,ref_db_I=None):
        '''Make a new SeqFeature'''

        # 1. Define the start/stop locations
        my_start_pos = SeqFeature.ExactPosition(start_pos_I)
        my_end_pos = SeqFeature.ExactPosition(stop_pos_I)

        # 2. Use the locations do define a FeatureLocation
        my_feature_location = FeatureLocation(my_start_pos,my_end_pos)

        # 4. Create a SeqFeature
        my_feature = SeqFeature.SeqFeature(location=my_feature_location,strand=strand_I,type=type_I,location_operator=location_operator_I,
                                           qualifiers=qualifiers_I,sub_features=sub_features_I,ref=ref_I,ref_db=ref_db_I)

        return my_feature;
        
    def _mutate_peptideFromMutationData(self,mutation_I,annotation_I=None,sequence_I=None,translation_table_I='Bacterial'):
        '''mutation dna, rna, and peptide sequences if the position given the reference genome
        input:
        mutation_I = {} mutation data, 
                        dictionary will be populated with different fields depending upon the mutation type
        annotation = genbank record for the annotation [SeqRecord]
        sequence = fasta record for the transcription/translation [SeqRecord]
        '''
        if not annotation_I:
            annotation_I = self.annotation;
        if not sequence_I:
            sequence_I = self.sequence;

        #parse the mutation dict
        if 'position' in mutation_I: mutation_position_I = mutation_I['position'];
        if 'type' in mutation_I: mutation_type_I = mutation_I['type'];
        if 'new_seq' in mutation_I: new_seq_I = mutation_I['new_seq'].lower();
        if 'size' in mutation_I: size_I = mutation_I['size'];
        elif 'new_seq' in mutation_I: size_I = len(new_seq_I);
        else: size_I = 0;
        if 'new_copy_number' in mutation_I: new_copy_number_I = mutation_I['new_copy_number'];

        snp_records = {};
        snp_records['gene'] = []
        snp_records['db_xref'] = []
        snp_records['locus_tag'] = []
        snp_records['EC_number'] = []
        snp_records['product'] = []
        snp_records['location'] = []
        snp_records['mutation_data'] = mutation_I;
        snp_records['dna_sequence_ref'] = ''
        snp_records['dna_sequence_new'] = ''
        snp_records['rna_sequence_ref'] = ''
        snp_records['rna_sequence_new'] = ''
        snp_records['peptide_sequence_ref'] = ''
        snp_records['peptide_sequence_new'] = ''
        snp_records['peptide_feature_position'] = None;
        snp_records['peptide_feature_ref'] = None;
        snp_records['peptide_feature_new'] = None;
        snp_records['peptide_feature_position'] = None;
        snp_records['peptide_feature_ref'] = None;
        snp_records['peptide_feature_new'] = None;
        snp_records['mutation_class'] = [];
        # find features in the coding region of the genome that bracket the mutation position
        for feature_cnt,feature in enumerate(annotation_I.features):
            if mutation_position_I in feature and feature.type == 'CDS':
                # extract out the original dna, rna, and peptide sequences
                dna,rna,peptide,has_stop_codon = self.transcribeAndTranslate_feature(sequence_I.seq,feature,table_I=translation_table_I);
                # mutate the sequence and adjust the feature size
                sequence_new,feature_new,mutation_new = self.mutate_sequence(sequence_I.seq,mutation_I,feature);
                if sequence_new is None:
                    break;
                # extract out the new dna, rna, and peptide sequences
                dna_new,rna_new,peptide_new,has_stop_codon_new = self.transcribeAndTranslate_feature(sequence_new,feature_new,table_I=translation_table_I);
                # convert to string output
                snp_records['gene'] = feature.qualifiers.get('gene')
                snp_records['location'] = ['coding'];
                snp_records['mutation_data'] = mutation_new;
                snp_records['dna_sequence_ref'] = self.convert_Seq2String(dna);
                snp_records['dna_sequence_new'] = self.convert_Seq2String(dna_new);
                snp_records['rna_sequence_ref'] = self.convert_Seq2String(rna);
                snp_records['rna_sequence_new'] = self.convert_Seq2String(rna_new);
                snp_records['peptide_sequence_ref'] = self.convert_Seq2String(peptide);
                snp_records['peptide_sequence_new'] = self.convert_Seq2String(peptide_new);
                # determine the mutation class
                mutation_class = {};
                mutation_class = self.classify_mutation(mutation_new,dna,rna,peptide,has_stop_codon,
                                                        dna_new,rna_new,peptide_new,has_stop_codon_new);
                #mutation_class = self.classify_mutation(mutation_I,dna,rna,peptide,has_stop_codon,
                #                                        dna_new,rna_new,peptide_new,has_stop_codon_new);
                snp_records.update(mutation_class);
                break;
        return snp_records;
    def compare2sequences(self,sequence_1,sequence_2):
        """Compare two sequences
        INPUT:
        sequence_1 = sequence object, list, or string
        sequence_2 = sequence object, list, or string
        OUTPUT:
        changed_pos = list of changed positions
        changed_sequences = list of changed sequences
        """
        changed_pos = [];
        old_sequences = [];
        new_sequences = [];
        cnt = None; # initialize the counter
        if len(sequence_1)>len(sequence_2):
            for cnt,p in enumerate(sequence_2):
                if p!=sequence_1[cnt]:
                    changed_pos.append(cnt+1);
                    old_sequences.append(sequence_1[cnt]);
                    new_sequences.append(p);
                    break;
                # check that a changed position was found
                # if not, the protein was truncated at the last position
            if not changed_pos and not cnt is None:
                changed_pos.append(cnt);
                old_sequences.append(sequence_1[cnt]);
                new_sequences.append('*');
        elif len(sequence_2)>len(sequence_1):
            for cnt,p in enumerate(sequence_1):
                if p!=sequence_2[cnt]:
                    changed_pos.append(cnt+1);
                    old_sequences.append(p);
                    new_sequences.append(sequence_2[cnt]);
                    break;
                # check that a changed position was found
                # if not, the protein was truncated at the last position
            if not changed_pos and not cnt is None:
                changed_pos.append(cnt);
                old_sequences.append(sequence_2[cnt]);
                new_sequences.append('*');
        else:
            for cnt,p in enumerate(sequence_1):
                if p!=sequence_2[cnt]:
                    changed_pos.append(cnt+1);
                    old_sequences.append(p);
                    new_sequences.append(sequence_2[cnt]);
        return changed_pos, old_sequences, new_sequences;

    def convert_ntSize2AASize(self,nt_size_I):
        '''convert nucleotide size to AA size'''
        aa_size = nt_size_I/3;
        aa_size_O = int(aa_size) + 1;
        return aa_size_O;

    def convert_genomePosition2FeaturePosition(self,mutation_position_I,feature_I):
        """convert genome position to feature position
        INPUT:
        mutation_position_I = Int, mutation position in the genome
        feature_I = SeqFeature object
        OUTPUT:
        feature_position_O = Int, mutation position in the feature
        """
        feature_position_O = None;
        feature_position_O = mutation_position_I - int(feature_I.location.start);
        return feature_position_O;

    def classify_mutation(self,mutation_I,
                          dna,rna,peptide,has_stop_codon,
                          dna_new,rna_new,peptide_new,has_stop_codon_new):
        '''classify the type of mutation
        Classifiers:
        missense = point mutation (nonsynonymous mutation that results in a different AA)
        synonymous = point mutation that does not result in a different AA
        nonsynonymous = 
        frameshift = a deletion that is not divisible by 3 that changes the downstream reading frame
        non-frameshift = a deletion that is divisible by 3 that does not change the downstream reading frame
        no stop codon = no stop codon found in the peptide
        nonsense = mutation that truncates the peptide at the site of the mutation
        truncated peptide = peptide with a shorter length than what would be expected
        '''

        #parse the mutation dict
        if 'position' in mutation_I: mutation_position_I = mutation_I['position'];
        if 'type' in mutation_I: mutation_type_I = mutation_I['type'];
        if 'new_seq' in mutation_I: new_seq_I = mutation_I['new_seq'].lower();
        if 'size' in mutation_I: size_I = mutation_I['size'];
        elif 'new_seq' in mutation_I: size_I = len(new_seq_I);
        else: size_I = 1;
        if 'new_copy_number' in mutation_I: new_copy_number_I = mutation_I['new_copy_number'];

        aa_size_I = self.convert_ntSize2AASize(size_I);

        mutation_class = {};
        mutation_class['dna_feature_position'] = None;
        mutation_class['dna_feature_ref'] = None;
        mutation_class['dna_feature_new'] = None;
        mutation_class['rna_feature_position'] = None;
        mutation_class['rna_feature_ref'] = None;
        mutation_class['rna_feature_new'] = None;
        mutation_class['peptide_feature_position'] = None;
        mutation_class['peptide_feature_ref'] = None;
        mutation_class['peptide_feature_new'] = None;
        mutation_class['mutation_class'] = [];
        if mutation_type_I=='SNP':
            pos,old,new=[],[],[];
            pos,old,new = self.compare2sequences(dna,dna_new);
            if not pos or not old or not new: return mutation_class;
            mutation_class['dna_feature_position'] = pos[0];
            mutation_class['dna_feature_ref'] = old[0];
            mutation_class['dna_feature_new'] = new_seq_I;
            pos,old,new=[],[],[];
            pos,old,new = self.compare2sequences(rna,rna_new);
            mutation_class['rna_feature_position'] = pos[0];
            mutation_class['rna_feature_ref'] = old[0];
            mutation_class['rna_feature_new'] = new[0];
            if pos and len(rna)>len(rna_new):
                mutation_class['mutation_class'].append('truncated transcript');
            pos,old,new=[],[],[];
            pos,old,new = self.compare2sequences(peptide,peptide_new);
            if pos and len(peptide)==len(peptide_new):
                mutation_class['mutation_class'].append('missense');
            elif pos and len(peptide_new)<pos[0]+aa_size_I:
                mutation_class['mutation_class'].append('nonsense');
            elif pos and len(peptide)!=len(peptide_new) and len(peptide)<len(peptide_new) and not has_stop_codon_new:
                mutation_class['mutation_class'].append('nonsynonymous');
                mutation_class['mutation_class'].append('no stop codon');
                mutation_class['mutation_class'].append('truncated peptide');
            elif pos and len(peptide)!=len(peptide_new) and len(peptide)<len(peptide_new):
                mutation_class['mutation_class'].append('nonsynonymous');
                mutation_class['mutation_class'].append('truncated peptide');
            elif pos and len(peptide)!=len(peptide_new) and not has_stop_codon_new:
                mutation_class['mutation_class'].append('nonsynonymous');
                mutation_class['mutation_class'].append('no stop codon');
                mutation_class['mutation_class'].append('change in peptide length');
            elif pos and len(peptide)!=len(peptide_new):
                mutation_class['mutation_class'].append('nonsynonymous');
                mutation_class['mutation_class'].append('change in peptide length');

            elif len(peptide)==len(peptide_new):
                mutation_class['mutation_class'].append('synonymous');
            elif len(peptide_new)<pos[0]+aa_size_I:
                mutation_class['mutation_class'].append('nonsense');
            elif len(peptide)!=len(peptide_new) and len(peptide)<len(peptide_new) and not has_stop_codon_new:
                mutation_class['mutation_class'].append('no stop codon');
                mutation_class['mutation_class'].append('truncated peptide');
            elif len(peptide)!=len(peptide_new) and len(peptide)<len(peptide_new):
                mutation_class['mutation_class'].append('truncated peptide');
            elif len(peptide)!=len(peptide_new) and not has_stop_codon_new:
                mutation_class['mutation_class'].append('no stop codon');
                mutation_class['mutation_class'].append('change in peptide length');
            elif len(peptide)!=len(peptide_new):
                mutation_class['mutation_class'].append('change in peptide length');

            else:
                mutation_class['mutation_class'].append('unclassified mutation');
            if pos:
                mutation_class['peptide_feature_position'] = pos[0];
                mutation_class['peptide_feature_ref'] = old[0];
                mutation_class['peptide_feature_new'] = new[0];
            else:
                mutation_class['peptide_feature_position'] = None;
                mutation_class['peptide_feature_ref'] = None;
                mutation_class['peptide_feature_new'] = None;
        elif mutation_type_I=='SUB':
            # check if the peptide sequence changed
            pos,old,new=[],[],[];
            pos,old,new = self.compare2sequences(dna,dna_new);
            if not pos or not old or not new: return mutation_class;
            mutation_class['dna_feature_position'] = pos[0];
            mutation_class['dna_feature_ref'] = old[0];
            mutation_class['dna_feature_new'] = new_seq_I;
            pos,old,new=[],[],[];
            pos,old,new = self.compare2sequences(rna,rna_new);
            mutation_class['rna_feature_position'] = pos[0];
            mutation_class['rna_feature_ref'] = old[0];
            mutation_class['rna_feature_new'] = new[0];
            if pos and len(rna)>len(rna_new):
                mutation_class['mutation_class'].append('truncated transcript');
            pos,old,new=[],[],[];
            pos,old,new = self.compare2sequences(peptide,peptide_new);
            if pos and len(peptide)==len(peptide_new):
                mutation_class['mutation_class'].append('missense');
            elif pos and len(peptide_new)<pos[0]+aa_size_I:
                mutation_class['mutation_class'].append('nonsense');
            elif pos and len(peptide)!=len(peptide_new) and len(peptide)<len(peptide_new) and not has_stop_codon_new:
                mutation_class['mutation_class'].append('nonsynonymous');
                mutation_class['mutation_class'].append('no stop codon');
                mutation_class['mutation_class'].append('truncated peptide');
            elif pos and len(peptide)!=len(peptide_new) and len(peptide)<len(peptide_new):
                mutation_class['mutation_class'].append('nonsynonymous');
                mutation_class['mutation_class'].append('truncated peptide');
            elif pos and len(peptide)!=len(peptide_new) and not has_stop_codon_new:
                mutation_class['mutation_class'].append('nonsynonymous');
                mutation_class['mutation_class'].append('no stop codon');
                mutation_class['mutation_class'].append('change in peptide length');
            elif pos and len(peptide)!=len(peptide_new):
                mutation_class['mutation_class'].append('nonsynonymous');
                mutation_class['mutation_class'].append('change in peptide length');
                
            elif len(peptide)==len(peptide_new):
                mutation_class['mutation_class'].append('synonymous');
            elif len(peptide_new)<pos[0]+aa_size_I:
                mutation_class['mutation_class'].append('nonsense');
            elif len(peptide)!=len(peptide_new) and len(peptide)<len(peptide_new) and not has_stop_codon_new:
                mutation_class['mutation_class'].append('no stop codon');
                mutation_class['mutation_class'].append('truncated peptide');
            elif len(peptide)!=len(peptide_new) and len(peptide)<len(peptide_new):
                mutation_class['mutation_class'].append('truncated peptide');
            elif len(peptide)!=len(peptide_new) and not has_stop_codon_new:
                mutation_class['mutation_class'].append('no stop codon');
                mutation_class['mutation_class'].append('change in peptide length');
            elif len(peptide)!=len(peptide_new):
                mutation_class['mutation_class'].append('change in peptide length');

            else:
                mutation_class['mutation_class'].append('unclassified mutation');
            if pos:
                mutation_class['peptide_feature_position'] = pos[0];
                mutation_class['peptide_feature_ref'] = old[0];
                mutation_class['peptide_feature_new'] = new[0];
            else:
                mutation_class['peptide_feature_position'] = None;
                mutation_class['peptide_feature_ref'] = None;
                mutation_class['peptide_feature_new'] = None;
        # check for frameshift mutations
        elif mutation_type_I=='DEL':
            pos,old,new=[],[],[];
            pos,old,new = self.compare2sequences(dna,dna_new);
            if not pos or not old or not new: return mutation_class;
            mutation_class['dna_feature_position'] = pos[0];
            mutation_class['dna_feature_ref'] = old[0];
            mutation_class['dna_feature_new'] = new[0];
            if size_I%3 == 0:
                mutation_class['mutation_class'].append('nonframeshift');
            else:
                mutation_class['mutation_class'].append('frameshift');
            pos,old,new=[],[],[];
            pos,old,new = self.compare2sequences(rna,rna_new);
            mutation_class['rna_feature_position'] = pos[0];
            mutation_class['rna_feature_ref'] = old[0];
            mutation_class['rna_feature_new'] = new[0];
            if pos and len(rna)-size_I>len(rna_new):
                mutation_class['mutation_class'].append('truncated transcript');
            pos,old,new=[],[],[];
            pos,old,new = self.compare2sequences(peptide,peptide_new);
            if pos and len(peptide)==len(peptide_new):
                mutation_class['mutation_class'].append('missense');
            elif pos and len(peptide_new)<pos[0]:
                mutation_class['mutation_class'].append('nonsense');
            elif pos and len(peptide)!=len(peptide_new) and len(peptide)-aa_size_I>len(peptide_new) and not has_stop_codon_new:
                mutation_class['mutation_class'].append('nonsynonymous');
                mutation_class['mutation_class'].append('no stop codon');
                mutation_class['mutation_class'].append('truncated peptide');
            elif pos and len(peptide)!=len(peptide_new) and len(peptide)-aa_size_I>len(peptide_new):
                mutation_class['mutation_class'].append('nonsynonymous');
                mutation_class['mutation_class'].append('truncated peptide');
            elif pos and len(peptide)!=len(peptide_new) and not has_stop_codon_new:
                mutation_class['mutation_class'].append('nonsynonymous');
                mutation_class['mutation_class'].append('no stop codon');
                mutation_class['mutation_class'].append('change in peptide length');
            elif pos and len(peptide)!=len(peptide_new):
                mutation_class['mutation_class'].append('nonsynonymous');
                mutation_class['mutation_class'].append('change in peptide length');
                
            elif len(peptide)==len(peptide_new):
                mutation_class['mutation_class'].append('synonymous');
            elif len(peptide_new)<pos[0]:
                mutation_class['mutation_class'].append('nonsense');
            elif len(peptide)!=len(peptide_new) and len(peptide)-aa_size_I>len(peptide_new) and not has_stop_codon_new:
                mutation_class['mutation_class'].append('no stop codon');
                mutation_class['mutation_class'].append('truncated peptide');
            elif len(peptide)!=len(peptide_new) and len(peptide)-aa_size_I>len(peptide_new):
                mutation_class['mutation_class'].append('truncated peptide');
            elif len(peptide)!=len(peptide_new) and not has_stop_codon_new:
                mutation_class['mutation_class'].append('no stop codon');
                mutation_class['mutation_class'].append('change in peptide length');
            elif len(peptide)!=len(peptide_new):
                mutation_class['mutation_class'].append('change in peptide length');

            else:
                mutation_class['mutation_class'].append('unclassified mutation');
            if pos:
                mutation_class['peptide_feature_position'] = pos[0];
                mutation_class['peptide_feature_ref'] = old[0];
                mutation_class['peptide_feature_new'] = new[0];
            else:
                mutation_class['peptide_feature_position'] = None;
                mutation_class['peptide_feature_ref'] = None;
                mutation_class['peptide_feature_new'] = None;
        elif mutation_type_I=='INS':
            pos,old,new=[],[],[];
            pos,old,new = self.compare2sequences(dna,dna_new);
            if not pos or not old or not new: return mutation_class;
            mutation_class['dna_feature_position'] = pos[0];
            mutation_class['dna_feature_ref'] = old[0];
            mutation_class['dna_feature_new'] = new[0];
            #new_seq = ''.join(new);
            #mutation_class['dna_feature_new'] = new_seq_I;
            if size_I%3 == 0:
                mutation_class['mutation_class'].append('nonframeshift');
            else:
                mutation_class['mutation_class'].append('frameshift');
            pos,old,new=[],[],[];
            pos,old,new = self.compare2sequences(rna,rna_new);
            mutation_class['rna_feature_position'] = pos[0];
            mutation_class['rna_feature_ref'] = old[0];
            mutation_class['rna_feature_new'] = new[0];
            if pos and len(rna)+size_I>len(rna_new)+1:
                mutation_class['mutation_class'].append('truncated transcript');
            pos,old,new=[],[],[];
            pos,old,new = self.compare2sequences(peptide,peptide_new);
            if pos and len(peptide)==len(peptide_new):
                mutation_class['mutation_class'].append('missense');
            elif pos and len(peptide_new)<pos[0]+aa_size_I:
                mutation_class['mutation_class'].append('nonsense');
            elif pos and len(peptide)!=len(peptide_new) and len(peptide)+aa_size_I>len(peptide_new) and not has_stop_codon_new:
                mutation_class['mutation_class'].append('nonsynonymous');
                mutation_class['mutation_class'].append('no stop codon');
                mutation_class['mutation_class'].append('truncated peptide');
            elif pos and len(peptide)!=len(peptide_new) and len(peptide)+aa_size_I>len(peptide_new):
                mutation_class['mutation_class'].append('nonsynonymous');
                mutation_class['mutation_class'].append('truncated peptide');
            elif pos and len(peptide)!=len(peptide_new) and not has_stop_codon_new:
                mutation_class['mutation_class'].append('nonsynonymous');
                mutation_class['mutation_class'].append('no stop codon');
                mutation_class['mutation_class'].append('change in peptide length');
            elif pos and len(peptide)!=len(peptide_new):
                mutation_class['mutation_class'].append('nonsynonymous');
                mutation_class['mutation_class'].append('change in peptide length');
                
            elif len(peptide)==len(peptide_new):
                mutation_class['mutation_class'].append('synonymous');
            elif len(peptide_new)<pos[0]+aa_size_I:
                mutation_class['mutation_class'].append('nonsense');
            elif len(peptide)!=len(peptide_new) and len(peptide)+aa_size_I>len(peptide_new) and not has_stop_codon_new:
                mutation_class['mutation_class'].append('no stop codon');
                mutation_class['mutation_class'].append('truncated peptide');
            elif len(peptide)!=len(peptide_new) and len(peptide)+aa_size_I>len(peptide_new):
                mutation_class['mutation_class'].append('truncated peptide');
            elif len(peptide)!=len(peptide_new) and not has_stop_codon_new:
                mutation_class['mutation_class'].append('no stop codon');
                mutation_class['mutation_class'].append('change in peptide length');
            elif len(peptide)!=len(peptide_new):
                mutation_class['mutation_class'].append('change in peptide length');

            else:
                mutation_class['mutation_class'].append('unclassified mutation');
            if pos:
                mutation_class['peptide_feature_position'] = pos[0];
                mutation_class['peptide_feature_ref'] = old[0];
                mutation_class['peptide_feature_new'] = new[0];
            else:
                mutation_class['peptide_feature_position'] = None;
                mutation_class['peptide_feature_ref'] = None;
                mutation_class['peptide_feature_new'] = None;
        elif mutation_type_I=='AMP':
            pos,old,new=[],[],[];
            pos,old,new = self.compare2sequences(dna,dna_new);
            if not pos or not old or not new: return mutation_class;
            mutation_class['dna_feature_position'] = pos[0];
            mutation_class['dna_feature_old'] = old[0];
            mutation_class['dna_feature_new'] = new[0];
            #amp_seq = self.convert_Seq2List(dna)[pos[0]-1:pos[0]-1+size_I];
            #new_seq_I = [];
            #for n in range(new_copy_number_I):
            #    new_seq_I.extend(amp_seq);
            # trim the sequence if it is too long
            #if len(new_seq_I)>100:
            #    mutation_class['dna_feature_new'] = new[0];
            #else:
            #    mutation_class['dna_feature_new'] = new_seq_I;
            new_seq_I = mutation_I['new_feature_sequence'];
            if len(new_seq_I)%3 == 0:
                mutation_class['mutation_class'].append('nonframeshift');
            else:
                mutation_class['mutation_class'].append('frameshift');
            pos,old,new=[],[],[];
            pos,old,new = self.compare2sequences(rna,rna_new);
            mutation_class['rna_feature_position'] = pos[0];
            mutation_class['rna_feature_ref'] = old[0];
            mutation_class['rna_feature_new'] = new[0];
            if pos and len(rna)+len(new_seq_I)>len(rna_new):
                mutation_class['mutation_class'].append('truncated transcript');
            aa_size_I = self.convert_ntSize2AASize(len(new_seq_I));
            pos,old,new=[],[],[];
            pos,old,new = self.compare2sequences(peptide,peptide_new);
            if pos and len(peptide)==len(peptide_new):
                mutation_class['mutation_class'].append('missense');
            elif pos and len(peptide_new)<pos[0]+aa_size_I:
                mutation_class['mutation_class'].append('nonsense');
            elif pos and len(peptide)!=len(peptide_new) and len(peptide)+aa_size_I>len(peptide_new) and not has_stop_codon_new:
                mutation_class['mutation_class'].append('nonsynonymous');
                mutation_class['mutation_class'].append('no stop codon');
                mutation_class['mutation_class'].append('truncated peptide');
            elif pos and len(peptide)!=len(peptide_new) and len(peptide)+aa_size_I>len(peptide_new):
                mutation_class['mutation_class'].append('nonsynonymous');
                mutation_class['mutation_class'].append('truncated peptide');
            elif pos and len(peptide)!=len(peptide_new) and not has_stop_codon_new:
                mutation_class['mutation_class'].append('nonsynonymous');
                mutation_class['mutation_class'].append('no stop codon');
                mutation_class['mutation_class'].append('change in peptide length');
            elif pos and len(peptide)!=len(peptide_new):
                mutation_class['mutation_class'].append('nonsynonymous');
                mutation_class['mutation_class'].append('change in peptide length');
                
            elif len(peptide)==len(peptide_new):
                mutation_class['mutation_class'].append('synonymous');
            elif len(peptide_new)<pos[0]+aa_size_I:
                mutation_class['mutation_class'].append('nonsense');
            elif len(peptide)!=len(peptide_new) and len(peptide)+aa_size_I>len(peptide_new) and not has_stop_codon_new:
                mutation_class['mutation_class'].append('no stop codon');
                mutation_class['mutation_class'].append('truncated peptide');
            elif len(peptide)!=len(peptide_new) and len(peptide)+aa_size_I>len(peptide_new):
                mutation_class['mutation_class'].append('truncated peptide');
            elif len(peptide)!=len(peptide_new) and not has_stop_codon_new:
                mutation_class['mutation_class'].append('no stop codon');
                mutation_class['mutation_class'].append('change in peptide length');
            elif len(peptide)!=len(peptide_new):
                mutation_class['mutation_class'].append('change in peptide length');

            else:
                mutation_class['mutation_class'].append('unclassified mutation');
            if pos:
                mutation_class['peptide_feature_position'] = pos[0];
                mutation_class['peptide_feature_ref'] = old[0];
                mutation_class['peptide_feature_new'] = new[0];
            else:
                mutation_class['peptide_feature_position'] = None;
                mutation_class['peptide_feature_ref'] = None;
                mutation_class['peptide_feature_new'] = None;
        elif mutation_type_I=='MOB':
            pos,old,new=[],[],[];
            pos,old,new = self.compare2sequences(dna,dna_new);
            if not pos or not old or not new: return mutation_class;
            mutation_class['dna_feature_position'] = pos[0];
            mutation_class['dna_feature_old'] = old[0];
            mutation_class['dna_feature_new'] = new[0];
            new_seq_I = mutation_I['new_feature_sequence'];
            if len(new_seq_I)%3 == 0:
                mutation_class['mutation_class'].append('nonframeshift');
            else:
                mutation_class['mutation_class'].append('frameshift');
            pos,old,new=[],[],[];
            pos,old,new = self.compare2sequences(rna,rna_new);
            mutation_class['rna_feature_position'] = pos[0];
            mutation_class['rna_feature_ref'] = old[0];
            mutation_class['rna_feature_new'] = new[0];
            if pos and len(rna)+len(new_seq_I)>len(rna_new):
                mutation_class['mutation_class'].append('truncated transcript');
            aa_size_I = self.convert_ntSize2AASize(len(new_seq_I));
            pos,old,new=[],[],[];
            pos,old,new = self.compare2sequences(peptide,peptide_new);
            if pos and len(peptide)==len(peptide_new):
                mutation_class['mutation_class'].append('missense');
            elif pos and len(peptide_new)<pos[0]+aa_size_I:
                mutation_class['mutation_class'].append('nonsense');
            elif pos and len(peptide)!=len(peptide_new) and len(peptide)+aa_size_I>len(peptide_new) and not has_stop_codon_new:
                mutation_class['mutation_class'].append('nonsynonymous');
                mutation_class['mutation_class'].append('no stop codon');
                mutation_class['mutation_class'].append('truncated peptide');
            elif pos and len(peptide)!=len(peptide_new) and len(peptide)+aa_size_I>len(peptide_new):
                mutation_class['mutation_class'].append('nonsynonymous');
                mutation_class['mutation_class'].append('truncated peptide');
            elif pos and len(peptide)!=len(peptide_new) and not has_stop_codon_new:
                mutation_class['mutation_class'].append('nonsynonymous');
                mutation_class['mutation_class'].append('no stop codon');
                mutation_class['mutation_class'].append('change in peptide length');
            elif pos and len(peptide)!=len(peptide_new):
                mutation_class['mutation_class'].append('nonsynonymous');
                mutation_class['mutation_class'].append('change in peptide length');
                
            elif len(peptide)==len(peptide_new):
                mutation_class['mutation_class'].append('synonymous');
            elif len(peptide_new)<pos[0]+aa_size_I:
                mutation_class['mutation_class'].append('nonsense');
            elif len(peptide)!=len(peptide_new) and len(peptide)+aa_size_I>len(peptide_new) and not has_stop_codon_new:
                mutation_class['mutation_class'].append('no stop codon');
                mutation_class['mutation_class'].append('truncated peptide');
            elif len(peptide)!=len(peptide_new) and len(peptide)+aa_size_I>len(peptide_new):
                mutation_class['mutation_class'].append('truncated peptide');
            elif len(peptide)!=len(peptide_new) and not has_stop_codon_new:
                mutation_class['mutation_class'].append('no stop codon');
                mutation_class['mutation_class'].append('change in peptide length');
            elif len(peptide)!=len(peptide_new):
                mutation_class['mutation_class'].append('change in peptide length');

            else:
                mutation_class['mutation_class'].append('unclassified mutation');
            if pos:
                mutation_class['peptide_feature_position'] = pos[0];
                mutation_class['peptide_feature_ref'] = old[0];
                mutation_class['peptide_feature_new'] = new[0];
            else:
                mutation_class['peptide_feature_position'] = None;
                mutation_class['peptide_feature_ref'] = None;
                mutation_class['peptide_feature_new'] = None;
        else:
            print('mutation type ' + mutation_type_I + ' not yet supported.');
        return mutation_class;

