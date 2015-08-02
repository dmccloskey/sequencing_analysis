from io_utilities.base_importData import base_importData
from io_utilities.base_exportData import base_exportData

class gd():
    def __init__(self,metadata_I,mutations_I,validation_I,evidence_I,mutationsAnnotated_I,mutationsFiltered_I):
        if metadata_I:
            self.metadata = metadata_I;
        else:
            metadata = [];
        if mutations_I:
            self.mutations=mutations_I;
        else:
            mutations = [];
        if validation_I:
            self.validation = validation_I;
        else:
            validation = [];
        if evidence_I:
            self.evidence = evidence_I;
        else:
            evidence = [];
        if mutationsAnnotated_I:
            self.mutationsAnnotated = mutationsAnnotated_I;
        else:
            mutationsAnnotated = [];
        if mutationsFiltered_I:
            self.mutationsFiltered = mutationsFiltered_I;
        else:
            mutationsFiltered = [];

    def import_gd(self,filename_I):
        """import and parse .gd file"""
        return
    def import_mutations(self,filename_I):
        """import mutations"""
        return
    def import_validation(self,filename_I):
        """import validation"""
        return
    def import_evidence(self,filename_I):
        """import validation"""
        return
    def import_mutationsAnnotated(self,filename_I):
        """import mutationsAnnotated"""
        return
    def import_mutationsFiltered(self,filename_I):
        """import mutationsFiltered"""
        return

    def export_metadata(self,filename_O):
        """export metadata"""
        return
    def export_mutations(self,filename_O):
        """export mutations"""
        return
    def export_validation(self,filename_O):
        """export validation"""
        return
    def export_evidence(self,filename_O):
        """export validation"""
        return
    def export_mutationsAnnotated(self,filename_O):
        """export mutationsAnnotated"""
        return
    def export_mutationsFiltered(self,filename_O):
        """export mutationsFiltered"""
        return

    def annotate_mutations(self):
        """annotate mutation from reference
        INPUT:
        """
        return
    def filter_mutations(self):
        """filter mutations based on user criteria
        INPUT:
        """