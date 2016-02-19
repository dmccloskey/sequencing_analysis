from io_utilities.base_importData import base_importData
from io_utilities.base_exportData import base_exportData
from sequencing_utilities import gdparse
from .genome_annotations import genome_annotations

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
    
    def format_mutationData(self,mutationData_I):
        """converts '{}' to {}"""
        for mutationData in mutationData_I:
            if 'mutation_data' in mutationData and type(mutationData['mutation_data'])==type('string'):
                mutationData['mutation_data'] = eval(mutationData['mutation_data']);
            if 'mutation_genes' in mutationData and type(mutationData['mutation_genes'])==type('string'):
                mutationData['mutation_genes'] = eval(mutationData['mutation_genes']);
            if 'mutation_locations' in mutationData and type(mutationData['mutation_locations'])==type('string'):
                mutationData['mutation_locations'] = eval(mutationData['mutation_locations']);
            if 'mutation_annotations' in mutationData and type(mutationData['mutation_annotations'])==type('string'):
                mutationData['mutation_annotations'] = eval(mutationData['mutation_annotations']);
            if 'mutation_frequency' in mutationData and type(mutationData['mutation_frequency'])==type('string'):
                mutationData['mutation_frequency'] = eval(mutationData['mutation_frequency']);
            if 'mutation_position' in mutationData and type(mutationData['mutation_position'])==type('string'):
                mutationData['mutation_position'] = eval(mutationData['mutation_position']);
        return mutationData_I;

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
        gd = gdparse.GDParser(file_handle=open(filename, 'r'))
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
        self.mutations = self.format_mutationData(io.data);
    def import_validation(self,filename_I):
        """import validation"""
        io = base_importData();
        io.read_csv(filename_I);
        self.validation = self.format_mutationData(io.data);
    def import_evidence(self,filename_I):
        """import evidence"""
        io = base_importData();
        io.read_csv(filename_I);
        self.evidence = self.format_mutationData(io.data);
    def import_mutationsAnnotated(self,filename_I):
        """import mutationsAnnotated"""
        io = base_importData();
        io.read_csv(filename_I);
        self.mutationsAnnotated = self.format_mutationData(io.data);
    def import_mutationsFiltered(self,filename_I):
        """import mutationsFiltered"""
        io = base_importData();
        io.read_csv(filename_I);
        self.mutationsFiltered = self.format_mutationData(io.data);

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

    def annotate_mutations(self,table_I='mutationsFiltered',ref_genome_I='U00096.2.gb',
                           ref_I = 'genbank',geneReference_I=None,biologicalmaterial_id_I='MG1655'):
        """annotate filtered mutation from reference
        INPUT:
        table_I = string name of the table to filter
        ref_genome_I = reference genome to use for the annotation
        ref_I = reference database
        geneReference_I = filename for the gene reference table
        biologicalmaterial_id_I = biologicalmatereial_id for the geneReference to use for the annotation (required to generate an ecocyc link)
        """
        genomeannotation = genome_annotations(ref_genome_I,ref_I,geneReference_I);

        # query mutation data:
        if table_I == 'mutationsFiltered':
            mutations = self.mutationsFiltered;
        elif table_I == 'mutationsLineage':
            mutations = self.mutationsLineage;
        elif table_I == 'mutationsEndpoints':
            mutations = self.mutationsEndpoints;
        else:
            mutations = self.mutations;
        mutation_data_O = [];
        for end_cnt,mutation in enumerate(mutations):
            data_tmp = {};
            # annotate each mutation based on the position
            annotation = {};
            annotation = genomeannotation._find_genesFromMutationPosition(mutation['mutation_data']['position']);
            data_tmp['mutation_genes'] = annotation['gene']
            data_tmp['mutation_locations'] = annotation['location']
            data_tmp['mutation_annotations'] = annotation['product']
            # generate a link to ecogene for the genes
            data_tmp['mutation_links'] = [];
            if genomeannotation.geneReference:
                for bnumber in annotation['locus_tag']:
                    if bnumber:
                        ecogenes = [];
                        ecogenes = genomeannotation._get_ecogenesByBiologicalmaterialIDAndOrderedLocusName(biologicalmaterial_id_I,bnumber);
                        if ecogenes:
                            ecogene = ecogenes[0];
                            ecogene_link = genomeannotation._generate_httplink2gene_ecogene(ecogene['ecogene_accession_number']);
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
                data_evidence_list = self._get_evidenceByPid(self.evidence,pid)
                # filter evidence based on user criteria
                data_evidence_filtered = [];
                data_evidence_filtered = self._filter_evidenceByPValueAndQualityAndFrequency(data_evidence_list,p_value_criteria=p_value_criteria,quality_criteria=quality_criteria,frequency_criteria=frequency_criteria);
                data_evidence_list.extend(data_evidence_filtered);
            if data_evidence_list: #check that filtered evidence was found
                data_O.append(data_mutations);
        #record the data
        self.mutationsFiltered = data_O;

    def _filter_mutationsByFrequency(self,data_I,frequency_criteria=0.1):
        """return data that is above the frequencing_criteria"""
        data_filtered = [];
        for d in data_I:
            if 'frequency' in d['mutation_data'] and d['mutation_data']['frequency'] >= frequency_criteria:
                data_filtered.append(d);
            #note: frequency is only provided for population resequences
            elif 'frequency' not in d['mutation_data']:
                data_filtered.append(d);
        return data_filtered;

    def _get_evidenceByPid(self,data_I,pid_I):
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


    def clear_data(self):
        """clear data"""
        del self.metadata[:];
        del self.evidence[:];
        del self.validation[:];
        del self.mutations[:];
        del self.mutationsFiltered[:];
        del self.mutationsAnnotated[:];

    def _make_mutationID(self,mutation_genes,mutation_type,mutation_position):
        '''return a unique mutation id string'''
        mutation_genes_str = '';
        for gene in mutation_genes:
            if mutation_genes_str is None: mutation_genes_str = 'unknown';
            mutation_genes_str = mutation_genes_str + gene + '-/-'
        mutation_genes_str = mutation_genes_str[:-3];
        mutation_id = mutation_type + '_' + mutation_genes_str + '_' + str(mutation_position);
        return mutation_id;

    def _make_mutationID2(self,mutation_genes,mutation_type,mutation_position,nt_ref,nt_new):
        '''return a unique mutation id string'''
        mutation_genes_str = '';
        for gene in mutation_genes:
            mutation_genes_str = mutation_genes_str + gene + '-/-'
        mutation_genes_str = mutation_genes_str[:-3];

        if nt_ref is None: nt_ref = '';
        if nt_new is None: nt_new = '';

        mutation_id = mutation_type + '_' + mutation_genes_str + '_' + nt_ref + str(mutation_position) + nt_new;
        return mutation_id;

    def _make_mutationGenesStr(self,mutation_genes):
        '''return a unique mutation genes list as a string'''
        mutation_genes_str = '';
        for gene in mutation_genes:
            #if mutation_genes_str is None: mutation_genes_str = 'unknown';
            mutation_genes_str = mutation_genes_str + gene + '-/-'
        mutation_genes_str = mutation_genes_str[:-3];
        return mutation_genes_str;

    def _make_mutationLocationsStr(self,mutation_locations):
        '''return a unique mutation Locations list as a string'''
        mutation_locations_str = '';
        for location in mutation_locations:
            #if mutation_locations_str is None: mutation_locations_str = 'unknown';
            mutation_locations_str = mutation_locations_str + location + '/'
        mutation_locations_str = mutation_locations_str[:-1];
        return mutation_locations_str;

    def _make_mutationClassStr(self,mutation_class):
        '''return a unique mutation Locations list as a string'''
        mutation_class_str = '';
        for location in mutation_class:
            #if mutation_class_str is None: mutation_class_str = 'unknown';
            mutation_class_str = mutation_class_str + location + '/'
        mutation_class_str = mutation_class_str[:-1];
        return mutation_class_str;