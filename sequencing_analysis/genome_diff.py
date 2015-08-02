from io_utilities.base_importData import base_importData
from io_utilities.base_exportData import base_exportData
from sequencing_utilities import gdparser
from Bio import SeqIO
from Bio import Entrez

class genome_diff():
    def __init__(self,metadata_I=None,mutations_I=None,validation_I=None,evidence_I=None,
                 mutationsAnnotated_I=None,mutationsFiltered_I=None,geneReference_I=None):
        if metadata_I:
            self.metadata = metadata_I;
        else:
            self.metadata = [];
        if mutations_I:
            self.mutations=mutations_I;
        else:
            self.mutations = [];
        if validation_I:
            self.validation = validation_I;
        else:
            self.validation = [];
        if evidence_I:
            self.evidence = evidence_I;
        else:
            self.evidence = [];
        if mutationsAnnotated_I:
            self.mutationsAnnotated = mutationsAnnotated_I;
        else:
            self.mutationsAnnotated = [];
        if mutationsFiltered_I:
            self.mutationsFiltered = mutationsFiltered_I;
        else:
            self.mutationsFiltered = [];

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
            self.geneReference = geneReference_I;
        else:
            self.geneReference = [];

    def import_biologicalMaterialGeneReferences(self, filename):
        '''import biological material gene references for annotations'''
        io = base_importData();
        io.read_csv(filename);
        self.geneReference=io.data;

    def import_gd(self, filename, experiment_id='', sample_name=''):
        """import and parse .gd file
        INPUT:
        filename = string, directory and filename of the .gd file

        OPTIONAL INPUT:
        the following are optional for analyzing a single sample,
        but required when analyzing multiple samples

        experiment_id = string, name of the experiment that generated the sample
        sample_name = string, name of the sample

        """
        gd = gdparse.GDParser(file_handle=open(filename, 'rb'))
        # extract out ids
        mutation_ids = [];
        mutation_ids = list(gd.data['mutation'].keys())
        parent_ids = [];
        for mid in mutation_ids:
            parents = [];
            parents = gd.data['mutation'][mid]['parent_ids'];
            parent_ids.extend(parents);
        # split into seperate data structures based on the destined table add
        metadata_data = [];
        mutation_data = [];
        evidence_data = [];
        validation_data = [];
        if gd.metadata:
            metadata_data.append({'experiment_id':experiment_id,
                                 'sample_name':sample_name,
                                 'genome_diff':gd.metadata['GENOME_DIFF'],
                                 'refseq':gd.metadata['REFSEQ'],
                                 'readseq':gd.metadata['READSEQ'],
                                 'author':gd.metadata['AUTHOR']});
        if gd.data['mutation']:
            for mid in mutation_ids:
                mutation_data.append({'experiment_id':experiment_id,
                                 'sample_name':sample_name,
                                 'mutation_id':mid,
                                 'parent_ids':gd.data['mutation'][mid]['parent_ids'],
                                 'mutation_data':gd.data['mutation'][mid]});
                                 #'mutation_data':json.dumps(gd.data['mutation'][mid])});
        if gd.data['evidence']:
            for pid in parent_ids:
                evidence_data.append({'experiment_id':experiment_id,
                                 'sample_name':sample_name,
                                 'parent_id':pid,
                                 'evidence_data':gd.data['evidence'][pid]});
                                 #'evidence_data':json.dumps(gd.data['evidence'][pid])});
        if gd.data['validation']:
            for mid in mutation_ids:
                validation_data.append({'experiment_id':experiment_id,
                                 'sample_name':sample_name,
                                 'validation_id':mid,
                                 'validation_data':gd.data['validation'][mid]});
                                 #'validation_data':json.dumps(gd.data['validation'][mid])});
        # add data to the database:
        self.metadata = metadata_data;
        self.mutations = mutation_data;
        self.evidence = evidence_data;
        self.validation = validation_data;
    def import_mutations(self,filename_I):
        """import mutations"""
        io = base_importData();
        io.read_csv(filename_I);
        self.mutations = io.data;
    def import_validation(self,filename_I):
        """import validation"""
        io = base_importData();
        io.read_csv(filename_I);
        self.validation = io.data;
    def import_evidence(self,filename_I):
        """import evidence"""
        io = base_importData();
        io.read_csv(filename_I);
        self.evidence = io.data;
    def import_mutationsAnnotated(self,filename_I):
        """import mutationsAnnotated"""
        io = base_importData();
        io.read_csv(filename_I);
        self.mutationsAnnotated = io.data;
    def import_mutationsFiltered(self,filename_I):
        """import mutationsFiltered"""
        io = base_importData();
        io.read_csv(filename_I);
        self.mutationsFiltered = io.data;

    def export_metadata(self,filename_O):
        """export metadata"""
        io = base_exportData(self.metadata);
        io.write_dict2csv(filename_O);
    def export_mutations(self,filename_O):
        """export mutations"""
        io = base_exportData(self.mutations);
        io.write_dict2csv(filename_O);
    def export_validation(self,filename_O):
        """export validation"""
        io = base_exportData(self.validation);
        io.write_dict2csv(filename_O);
    def export_evidence(self,filename_O):
        """export evidence"""
        io = base_exportData(self.evidence);
        io.write_dict2csv(filename_O);
    def export_mutationsAnnotated(self,filename_O):
        """export mutationsAnnotated"""
        io = base_exportData(self.mutationsAnnotated);
        io.write_dict2csv(filename_O);
    def export_mutationsFiltered(self,filename_O):
        """export mutationsFiltered"""
        io = base_exportData(self.mutationsFiltered);
        io.write_dict2csv(filename_O);

    def annotate_mutations(self,mutationsFiltered_I=True,ref_genome_I='U00096.2.gb',biologicalmaterial_id_I='MG1655'):
        """annotate filtered mutation from reference
        INPUT:
        mutationsFiltered_I = boolean, annotates the filtered mutations if true, annotates all mutations if false
        ref_genome_I = reference genome to use for the annotation
        biologicalmaterial_id_I = biologicalmatereial_id for the geneReference to use for the annotation (required to generate an ecocyc link)
        """
        record = SeqIO.read(ref_genome_I,'genbank');
        # query mutation data:
        if mutationsFiltered_I:
            mutations = self.mutationsFiltered;
        else:
            mutations = self.mutations;
        mutation_data_O = [];
        for end_cnt,mutation in enumerate(mutations):
            # annotate each mutation based on the position
            annotation = {};
            annotation = self._find_genesFromMutationPosition(mutation['mutation_data']['position'],record);
            data_tmp['mutation_genes'] = annotation['gene']
            data_tmp['mutation_locations'] = annotation['location']
            data_tmp['mutation_annotations'] = annotation['product']
            # generate a link to ecogene for the genes
            data_tmp['mutation_links'] = [];
            if self.geneReference:
                for bnumber in annotation['locus_tag']:
                    if bnumber:
                        ecogenes = [];
                        ecogenes = self._get_ecogenesByBiologicalmaterialIDAndOrderedLocusName('MG1655',bnumber);
                        if ecogenes:
                            ecogene = ecogenes[0];
                            ecogene_link = self._generate_httplink2gene_ecogene(ecogene['ecogene_accession_number']);
                            data_tmp['mutation_links'].append(ecogene_link)
                        else: print('no ecogene_accession_number found for ordered_locus_location ' + bnumber);
            data_tmp['experiment_id'] = mutation['experiment_id'];
            data_tmp['sample_name'] = mutation['sample_name'];
            frequency = 1.0;
            if 'frequency' in mutation['mutation_data']:
                frequency = mutation['mutation_data']['frequency'];
            data_tmp['mutation_frequency'] = frequency
            data_tmp['mutation_position'] = mutation['mutation_data']['position']
            data_tmp['mutation_type'] = mutation['mutation_data']['type']
            data_tmp['mutation_data'] = mutation['mutation_data'];
            mutation_data_O.append(data_tmp);
        #record the data
        self.mutationsAnnotated = mutation_data_O;

    def filter_mutations_population(self,p_value_criteria=0.01,quality_criteria=6.0,frequency_criteria=0.1):
        """filter mutations based on user criteria
        INPUT:
        p_value_criteria=0.01 (default)
        quality_criteria=6.0 (default)
        frequency_criteria=0.1 (default)
        """
        data_O = [];
        #query mutation data filtered by frequency
        data_mutations_list = [];
        data_mutations_list = self._filter_mutationsByFrequency(self.mutations,frequency_criteria);
        for data_mutations in data_mutations_list:
            data_evidence_list = [];
            for pid in data_mutations['parent_ids']:
                # get evidence for the pid
                data_evidence_list = [];
                data_evidence_list = self._get_evidenceByPid(self.evidence,pid_I)
                # filter evidence based on user criteria
                data_evidence_filtered = [];
                data_evidence_filtered = self._filter_evidenceByPValueAndQualityAndFrequency(data_evidence_list,p_value_criteria=p_value_criteria,quality_criteria=quality_criteria,frequency_criteria=frequency_criteria);
                data_evidence_list.extend(data_evidence_filtered);
            if data_evidence_list: #check that filtered evidence was found
                data_O.append(data_mutations);
        #record the data
        self.mutationsFiltered = data_O;

    def _filter_mutationsByFrequencing(self,data_I,frequencing_criteria=0.1):
        """return data that is above the frequencing_criteria"""
        data_filtered = [];
        for d in data_I:
            if 'frequency' in d['mutation_data'] and d['mutation_data']['frequency'] >= frequency_criteria:
                data_filtered.append(d);
            #note: frequency is only provided for population resequences
            elif 'frequency' not in d['mutation_data']:
                data_filtered.append(d);
        return data_filtered;

    def _get_evidenceByPid(data_I,pid_I):
        """return data with matching pid"""
        data_evidence_list = [];
        data_evidence_dict = {};
        for d in data_I:
            if d['parent_id']==pid_I:
                data_evidence_dict = d;
                data_evidence_list.append(d);
        return data_evidence_list;

    def _filter_evidenceByPValueAndQualityAndFrequency(self,data_I,p_value_criteria=0.01,quality_criteria=6.0,frequency_criteria=0.1):
        """return data this is above the p-value, quality, and frequency criteria"""
        
        data_filtered = [];
        data_filtered_dict = {};
        for d in data_I:
            if d['evidence_data']['type'] == 'RA':
                #population only
                if 'quality' in d['evidence_data'] and \
                    'bias_p_value' in d['evidence_data'] and \
                    'fisher_strand_p_value' in d['evidence_data'] and \
                    'frequency' in d['evidence_data'] and \
                    d['evidence_data']['frequency'] >= frequency_criteria and \
                    d['evidence_data']['quality'] >= quality_criteria and \
                    d['evidence_data']['bias_p_value'] <= p_value_criteria and \
                    d['evidence_data']['fisher_strand_p_value'] <= p_value_criteria:
                    data_filtered_dict = d;
                    data_filtered.append(d);
                #population and isolate
                elif 'quality' in d['evidence_data'] and \
                    'frequency' in d['evidence_data'] and \
                    d['evidence_data']['frequency'] >= frequency_criteria and \
                    d['evidence_data']['quality'] >= quality_criteria:
                    data_filtered_dict = d;
                    data_filtered.append(d);
            elif d['evidence_data']['type'] == 'JC':
                data_filtered_dict = d;
                data_filtered.append(d);
            elif d['evidence_data']['type'] == 'MC':
                data_filtered_dict = d;
                data_filtered.append(d);
            elif d['evidence_data']['type'] == 'UN':
                data_filtered_dict = d;
                data_filtered.append(d);
            else:
                print('mutation evidence of type ' + d['evidence_data']['type'] +\
                    ' has not yet been included in the filter criteria');
        return data_filtered;

    def _find_genesFromMutationPosition(self,mutation_position_I,record_I):
        '''find genes at the position or closest to the position given the reference genome'''
        #input:
        # mutation_position_I = mutation position [int]
        # record = genbank record [SeqRecord]
        snp_records = {};
        snp_records['gene'] = []
        snp_records['db_xref'] = []
        snp_records['locus_tag'] = []
        snp_records['EC_number'] = []
        snp_records['product'] = []
        snp_records['location'] = []
        # find features in the coding region of the genome that bracket the mutation position
        for feature_cnt,feature in enumerate(record_I.features):
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
            for feature_cnt,feature in enumerate(record_I.features):
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
            for feature_cnt,feature in enumerate(record_I.features):
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
    def _find_genesInRegion(self,start_I,stop_I,record_I):
        '''find genes in the start and stop region of the genome
        INPUT:
        mutation_position_I = mutation position [int]
        record_I = genbank record [SeqRecord]
        '''
        data_O = [];
        #extract all features within the start and stop region
        features = [f for f in record_I.features if start_I <= f.location.start.position and stop_I <= f.location.end.position]
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
    def _get_ecogenesByBiologicalmaterialIDAndOrderedLocusName(biologicalmaterial_id_I,ordered_locus_name_I):
        """return the ecogenes that match the biologicalmaterial_id and ordered_locus_name"""
        ecogenes = [];
        for d in self.geneReferences:
            if d['biologicalmaterial_id'] ==  biologicalmaterial_id_I and d['ordered_locus_name'] == ordered_locus_name_I:
                ecogenes.append(d);
        return ecogenes;

    def clear_data(self):
        """clear data"""
        del self.metadata[:];
        del self.evidence[:];
        del self.validation[:];
        del self.mutations[:];
        del self.mutationsFiltered[:];
        del self.mutationsAnnotated[:];