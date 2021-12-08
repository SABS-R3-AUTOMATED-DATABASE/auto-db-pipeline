'''

Database_interface.py

A module to provide a python API for the Antibody Database.


@author: James Dunbar
@contact: james.dunbar-[at]-dtc.ox.ac.uk

'''

# python
import os, sys

# ABDB
from ABDB.AbPDB.AntibodyParser import AntibodyParser
from ABDB.Annotate.annotate import get_alignment_dict
from ABDB.AB_Utils.AB_Utils import uniq, aa1, tuple_interpret, sci_to_common_names
from ABDB.AB_Utils.region_definitions import annotate_regions
from ABDB.AB_Utils.canonicals import assign_canonical
from ABDB import database_path, abysis, numbering_software, allow_online, anarci_available

class Database(object):
    """
    The database class.
    
    This is the back-end parser and interface for SAbDab.
     
    The "fetch" method is used to get an object that contains the details about a structure.
    
    e.g. API useage
    >>> from ABDB import database # database is initialised when ABDB is called.
    >>> p = database.fetch("12e8")
    >>> p.get_species()
    'MUS MUSCULUS'
    
    """
    def __init__(self, dbpath=""):
        if dbpath:
            self.dbpath = dbpath
        else:
            print("Path to database not provided. Run the setup script.", file=sys.stderr) 
            raise
              
        if os.path.exists(self.dbpath):
            pass
        else:
            print("Database path %s was not found." % self.dbpath, file=sys.stderr)
            self.loaded = False
            return

        # These are the fields we expect in the database summary.
        # Some are fab specific, some are for the whole pdb file.

        # Fields should be added here if they appear the the summary file.
        self.pdbfields = ["pdb", "short_header", "date", "compound", "organism", "authors", "resolution", "method", "r_free", "r_factor", "affinity", "temperature","affinity_method", "delta_g", "pmid"] 
    
        self.fabfields = ["Hchain", "Lchain", "model", "antigen_chain", "antigen_type", "antigen_het_name", "antigen_name", "scfv", "engineered", "heavy_subclass", "light_subclass", "light_ctype", "heavy_species", "light_species", "antigen_species"]

        self.banned_list = [ '4hjg', '4hjj', '4org', '4od1', '4ocw', '4ocs', '4od3', '4odh', '4ocr', '3upc', '4unt', '3eot', '1sjv', '1h8o', '1h8s', '1fn4', '3mlt', '1ngx', '3mls', '3mlr', '3mlu','3mlv', '3dvg', '3dvn', '1oau', '2gki', '1lmk', '1nqb', '1moe', '4k3e', '4k3d', '4kro','5dt1', '4dqo', '4hqq','4isv','4jm4','4lrt','4lrn','4lsp','4lsq','4lsr','4nuj','4ob5','4osu','4pub','4x8j','4hx2','4yaq','3c08','3d0l','3drq','3lrs','3mug','3na9','3qeg','3qhz','3u1s','3u2s','3u4e','1fgv','1mie','1mj7','1rzf','1rzg','2a77','2dtm','2hff','2iqa','5e7b','5ihu','5e99','5ijv','5ilt','3bj9','5hcg','5fhx','5fcs','5gry','5gru','5drn','1efq','7fab']

        self.fields = set(self.pdbfields + self.fabfields)

        self.db_summary = {}
        self.no_fvs = 0

        # read the summary file - basic information about each fab structure.
        if self._read_db_summary():
            self.loaded = False
        else:
            self.loaded = True

        # attempt to read the in-house summary file - basic information about each in-house fab structure.
        if self._read_db_summary(inhouse=True):
            self.loaded_inhouse = False
        else:
            self.loaded_inhouse = True

    def __iter__(self):
        for pdb in self.db_summary:
            yield pdb
            
    def __str__(self):
        return """
        SAbDab                                         \\\    //
        The OPIG Structural Antibody Database           \\\  //
        Authors: James Dunbar and Konrad Krawczyk 2013-15 ||
        Contributors: Jinwoo Leem                         ||
        Supervisor: Charlotte Deane
        
        Contact: Charlotte Deane (deane@stats.ox.ac.uk)
        
        Please cite:
        Dunbar J, Krawczyk K et al (2014) Nucleic Acids Res. 
        
        Last update to database found at %s:
        %s
        
        This database contains %d structures from the PDB containing %d fv regions.
        """ % (self.dbpath, self.get_update_date(), len(self.db_summary), self.no_fvs)


    # public methods
    def fetch(self, pdb):
        """
        Use this function to get the details for a structure in the database.
        
        If the code is in our database, a L{pdb_details<pdb_details>} object will be returned.
        
        @param pdb: A four alphanumeric pdb code.
        @type pdb: C{str}
        """

        pdb = pdb.lower()

        if pdb in self:
            details = self.db_summary[pdb]

            # Create a details object
            pdb_object = PDB_details(pdb,database=self)

            #####################################
            # Details relating to the whole pdb #
            #####################################
            # read from details dictionary
            pdb_object._set_short_header(details["short_header"])
            pdb_object._set_compound(details["compound"])
            pdb_object._set_date(details["date"])
            pdb_object._set_authors(details["authors"])
            pdb_object._set_organism(details["organism"])
            pdb_object._set_method(details["method"])
            pdb_object._set_resolution(details["resolution"])
            pdb_object._set_r_factor(details["r_factor"])  # added 190413
            pdb_object._set_r_free(details["r_free"])  # added 190413
            pdb_object._set_affinity(details["affinity"])  # added 220513
            pdb_object._set_delta_g(details["delta_g"])  # added 220513
            pdb_object._set_pmid(details["pmid"])  # added 220513
            pdb_object._set_affinity_method(details["affinity_method"])  # added 300913
            pdb_object._set_temperature(details["temperature"])  # added 300913
            try:
                pdb_object._set_inhouse(details["inhouse"])  # added 280515. Work around with "try" required for pickled db
            except KeyError:
                pdb_object._set_inhouse(False)
              
            # set paths to larger data in the database (Sequences, numbering, structure etc) 
            # check existence of files but do not fetch until required.
            # This makes loadall feasible.  
            pdb_object._set_headerpath(self._get_headerpath(pdb))
            pdb_object._set_filepath(self._get_filepath(pdb))                        
            pdb_object._set_annotationpath(self._get_annotationpath(pdb))
            pdb_object._set_sequencepath(self._get_sequencepath(pdb))
            pdb_object._set_imgtpath(self._get_imgtpath(pdb))
            pdb_object._set_abanglepath(self._get_abanglepath(pdb))
            pdb_object._set_summarypath(self._get_summarypath(pdb))
            pdb_object._set_clusteringpath(self._get_clusteringpath(pdb))
            
            #######################################
            # Details relating to specific chains #
            #######################################
            # fabs             
            pdb_object._set_fabs(details["fabs"])  # set fab chains and their details - done in the pdb_details class. 

            # antigens
            pdb_object._set_antigens([fab.antigen for fab in pdb_object.get_fabs() if fab.antigen != None])  # set antigen chains and their details - done in the pdb_details class.

            pdb_object._set_engineered()  # inherits from its fabs.

            return pdb_object

        else:
#            print "pdb %s not found in our database" % pdb
            return None

    def get_summary(self):
        """
        Get a summary of the entries in the database.

        @return: The database summary dictionary. Keys are the pdbs. Values are pdb details and a further dictionary for the fabs contained within the pdb file.  
        @rtype: C{dict}
        """
        return self.db_summary


    def get_update_date(self):
        """
        Get the date of the last update to the loaded database.
        @return: The date of the last update to the database.
        @rtype: C{str}
        """
        if hasattr(self, "update_date"): 
            return self.update_date
        try:
            logs = os.listdir(os.path.join(self.dbpath, "Log"))
            self.update_date = " ".join(sorted([log.split("_") for log in logs] , key=lambda l: (int(l[3]), "JANFEBMARAPRMAYJUNJULAUGSEPOCTNOVDEC".index(l[2].upper()), int(l[1]), int(l[4][:8].replace(":", ""))))[-1]).split(".")[1]
        except OSError:
            self.update_date = "unknown"
        finally:
            return self.update_date
    
    def get_latest_structures(self):
        """
        Get a list of the latest pdb structures to be added to the database.
        """
        try:
            latest=[]
            with open(os.path.join(self.dbpath, "Log", "dbupdate."+"_".join(self.get_update_date().split())+".log") ) as logfile:
                for line in logfile:
                    if "updated" in line:
                        pdb=line.split()[0]
                        latest.append(self.fetch(pdb))
            return latest
        except IOError:
            print("Could not open logfile", file=sys.stderr)
            return []
            
    def deposition_plot(self, png=""):
        """
        Generate a plot of the number of antibody structures in the PDB by year. 
        This uses the the deposition date from PDB files.
        
        A matplotlib plot is produced which can then be saved.
        
        Note that the numbers per year seem to decline in the current and previous year.
        My guess is that this is because the PDB may hold back structures until checks are made - therefore the latest update may have deposition dates from the previous year (or two)

        @requires: matplotlib, numpy
        """
        try:
            import numpy as np
            import matplotlib
            if png:matplotlib.use('Agg')
            import matplotlib.pyplot as plt
        except ImportError:
            print("numpy and matplotlib required to show deposition plot. Please install.", file=sys.stderr)
            return
        
        # apologies for the horrible line. Will work upto year 2050...
        dates = sorted([ (int(l[0]), int(l[1]), int("19" + l[2])) if int(l[2]) > 50 else (int(l[0]), int(l[1]), int("20" + l[2])) for l in [ self.db_summary[pdb]["date"].split("/") for pdb in self if len(pdb)==4] ] , key=lambda x: (x[2], x[0], x[1]) )
        year_bins = {}
        for date in dates:
            try:
                year_bins[ date[2] ] += 1
            except KeyError:
                year_bins[ date[2] ] = 1

        frequency_year = []
        for year in range(dates[0][-1], dates[-1][-1] + 1):
            if year in year_bins:
                frequency_year.append((year, year_bins[year]))
            else:
                frequency_year.append((year, 0)) 
        
        cum_sum, n = [], 0
        for m in frequency_year:
            n += m[1]
            cum_sum.append((m[0], n)) 

        ind = np.arange(len(cum_sum))
        width = 0.35
        fig = plt.figure(figsize=(10, 6.5))
        ax = fig.add_subplot(111)
        rects = ax.bar(ind, [m[1] for m in cum_sum], width, color='r')
        lines = ax.plot(ind[:-1], [y[1] for y in frequency_year][:-1], 'bs-', ms=5)
        

        ax.set_ylabel('Number of Antibody Structures')
        ax.set_title('Number of Antibody Structures in the PDB by year\n(as of %s)' % self.get_update_date())
        
        xticks = ax.set_xticks(ind + width)
        xtick_labels = ax.set_xticklabels([str(m[0]) for m in cum_sum] , rotation='vertical')
        ax.legend((rects[0], lines[0]), ("Cumulative Sum", "Depositions per year"), loc='best')
        if png:
            fig.savefig(png)
        else:
            fig.show()
    
    def get_template(self, seq, region=[], strict=False, fallback=False, n=10, use_similarity=False):
        """
        Get the top n templates by sequence identity for your sequence.
        
        @param seq: The sequence you wish to model. Use the chain break symbol "/" to submit heavy and light chains together. e.g QEIDLF/SDFMMSD
        @param region: The region over which you wish to calculate sequence identity.
        @param n: The number of templates that should be returned. Default is 10
        @param strict: If True only templates that have 100% coverage over the region specified (or variable region) will be returned.
        @param fallback: If abnum fails to number the sequence try to do so with muscle number.
        @param use_similarity: If True, use sequence similarity instead of sequence identity. Sequence similarity calculated using blosum62
        
        
        @return: A list of templates [ ((pdb, chain, model), identity), ....]  or for h/l [ ((pdb, H, L, model), identity), ....]  AND the numbered sequence that was submitted.

        ["vl","vh","framework","interface","cdrs","cdrh1","cdrh2","cdrh3","cdrl1","cdrl2","cdrl3"]
        """    

        from ABDB.Annotate import online, abnum, muscle_number, anarci  # import numbering methods 
       
        if anarci_available:
            number = anarci
        elif abysis:
            number = abnum
        elif numbering_software:
            raise Exception("Unimplemented")
        elif allow_online:
            number = online
            
        
        sequence_dict = {}
        submitted_seqs = seq.split("/")
        for s in submitted_seqs:
            self._validate_sequence(s)
            if anarci_available:
                numbering, chain_type = number(s.upper(), scheme=self.get_numbering_scheme()) # Use the default numbering scheme for the database.
            else:
                numbering, chain_type = number(s.upper())
            if not numbering and fallback:  # abnum has failed but we will try to use muscle_number
                numbering, chain_type = muscle_number(s.upper())
            if not chain_type:
                raise AssertionError("Unable to number the given sequence")
            if chain_type in sequence_dict:
                raise AssertionError("Two %s chains submitted." % chain_type)
            if len(submitted_seqs) == 2 and sequence_dict and chain_type != "L":
                raise AssertionError("Light sequence is not an antibody light chain.")
            if len(submitted_seqs) == 2 and not sequence_dict and chain_type != "H":
                raise AssertionError("Heavy sequence is not an antibody heavy chain.")
            sequence_dict[chain_type] = dict(numbering)
         
        identities = self.get_identities(sequence_dict, region=region, strict=strict, use_similarity=use_similarity)

        numbered_sequence = {} 
        for chain in sequence_dict:
            numbered_sequence[chain] = sorted(list(sequence_dict[chain].items()), key=lambda x:x[0])

        # only return templates that are > 0% identical (e.g if you have defined a scrict cdr)
        top_temps = [ temp  for temp in sorted(list(identities.items()), key=lambda x: x[1] , reverse=True)[:n] if temp[1] > 0]
                            
        return top_temps, numbered_sequence

    # JD 13-03-15
    def get_identities(self,sequence_dict, region=[], strict=False, use_similarity=False):
        """
        Get the identities of a numbered sequence to the variable regions in the database.

        @param sequence_dict: A dictionary containing the numbered sequence position-aa dictionary. Chain type (H or L) as keys.
        @param region: The region over which you wish to calculate sequence identity.
        @param strict: If True only templates that have 100% coverage over the region will be given a sequence identity greater than 0.
        @param use_similarity: If True, use sequence similarity instead of sequence identity. Sequence similarity calculated using blosum62

        @return A dictionary containing the identities to the database entries.
        """
        from ABDB.AB_Utils.calculations import identity, similarity
        # Load the sequence database.
        # This will the default numbering scheme for the database. 
        # This can be changed by calling set_numbering_scheme *before* the first call of this method.
        self._get_numbered_sequences() 
        if not self.sequences:
            raise AssertionError("Sequences did not load - failed to find template")
            return

        # Check the regions.
        regions = []
        if region:
            for r in region:
                assert r.lower() in ["vh", "vl", "framework", "cdrs", "cdrh1", "cdrh2", "cdrh3", "cdrl1", "cdrl2", "cdrl3", "hframework", "lframework", "hcdrs", "lcdrs"], "Unrecognised/implemented region %s" % r
                regions.append(r.lower())
        else:
            for chain_type in sequence_dict:
                regions.append("v" + chain_type)

        assert len(sequence_dict) < 3, "Unexpected number of sequences identified"
        if len(list(sequence_dict.keys())) == 2:
            search_for = "HL"
        else:
            search_for = list(sequence_dict.keys())[0]
            
        if (set(["cdrh1", "cdrh2", "cdrh3"]) & set(regions)) and "H" not in search_for:
           raise AssertionError("Templates for heavy CDRs cannot be found without submitting a heavy sequence.") 
        if (set(["cdrl1", "cdrl2", "cdrl3"]) & set(regions)) and "L" not in search_for:
           raise AssertionError("Templates for light CDRs cannot be found without submitting a light sequence.") 

        if use_similarity:
            compare = similarity
        else:
            compare = identity

        identities = {}
        for pdb in self.sequences:
            for fab in self.sequences[pdb]:
                temp_seq = self.sequences[pdb][fab]
                if any(1 for s in search_for if s not in temp_seq): continue
                if search_for == "HL":
                    identities[(pdb, fab.id[1], fab.id[2], fab.id[3])] = max(compare(sequence_dict, temp_seq, region=regions, strict=strict, scheme=self.get_numbering_scheme(), definition=self.get_region_definition()), 0)
                elif search_for == "H":
                    identities[(pdb, fab.id[1], fab.id[3])] = max(compare(sequence_dict, temp_seq, region=regions, strict=strict, scheme=self.get_numbering_scheme(), definition=self.get_region_definition()), 0)
                elif search_for == "L": 
                    identities[(pdb, fab.id[2], fab.id[3])] = max(compare(sequence_dict, temp_seq, region=regions, strict=strict, scheme=self.get_numbering_scheme(), definition=self.get_region_definition()), 0)

        return identities

    # JD 13-03-15
    def get_identical(self, sequence_dict, region=[], threshold=0.99,strict=False, use_similarity=False, get_fabs=False ):
        """
        Get identical (or within a sequence identity threshold) structures to your sequence.

        @param sequence_dict: A dictionary containing the numbered sequence as a position-aa dictionary. Chain type (H or L) as keys.
        @param region: The region over which you wish to calculate sequence identity.
        @param threshold: The threshold at which to call things the same (default is >99%). Changes to blosum62 score if use_similarity is true
        @param strict: If True only templates that have 100% coverage over the region will be given a sequence identity greater than 0.
        @param use_similarity: If True, use sequence similarity instead of sequence identity. Sequence similarity calculated using blosum62
        @param get_fabs: Fetch these structures from the database? These fab objects become the values of the returned dictionary

        """
        if get_fabs: # Pre-check
            assert "H" in list(sequence_dict.keys()) and "L" in list(sequence_dict.keys()), "Get fabs can only be used with paired chains currently"

        # Get the identities to the database structures.
        identities = self.get_identities(sequence_dict, region=region, strict=strict, use_similarity=use_similarity)

        # Filter by the threshold and return. 
        if get_fabs:
            return dict( (identifier, self.fetch(identifier[0])[identifier[1]+identifier[2]] ) for identifier, identity in identities.items() if identity > threshold )
        else:
            return dict( (identifier, identity) for identifier, identity in identities.items() if identity > threshold )


    def set_numbering_scheme(self, scheme="chothia"):
        """
        Use the specified numbering scheme 
        """
        assert scheme.lower() in ["chothia", "kabat", "martin", "imgt", "wolfguy"], "Numbering scheme '%s' not recognised"%scheme
        self.numbering_scheme = scheme.lower()
    
    def get_numbering_scheme(self):
        """
        Get the numbering scheme used by the database.
        This can be changed using the set_numbering_scheme method
        """
        if not hasattr(self, "numbering_scheme"):
            self.set_numbering_scheme()
        elif not self.numbering_scheme:
            self.set_numbering_scheme()
        return self.numbering_scheme

    def set_region_definition(self, region_definition="chothia"):
        """
        Use the specified definition for regions (e.g. CDRs)
        """
        assert region_definition.lower() in ["chothia", "kabat", "north","contact","imgt", 'wolfguy'], "Definition name '%s' not recognised"%region_definition
        self.region_definition = region_definition.lower()
    
    def get_region_definition(self):
        """
        Get the region definition used by the database.
        This can be changed using the set_region_definition method
        """
        if not hasattr(self, "region_definition"):
            self.set_region_definition()
        elif not self.region_definition:
            self.set_region_definition()
        return self.region_definition
       
    # private methods
    def _read_db_summary(self, inhouse=False):
        """
        Read the summary file from the database.
        
        @change: Check to see whether the summary file found in the loaded database has all the fields expected by the interface.
        If not, the user is told and None is passed to the fields without information.
        """

        if inhouse:
            filepath = os.path.join(self.dbpath, "Summaries", "inhouse_db_summary.dat")
        else:
            filepath = os.path.join(self.dbpath, "Summaries", "db_summary.dat") 

        try:
            summary_file = open(filepath)
        except IOError:
            if not inhouse:
                print("No database summary file found", file=sys.stderr)
            return 1
      
        lines = summary_file.readlines()
        
        if not lines:
            print("No database summary found", file=sys.stderr)
            return 1
        
        self.no_fvs += len(lines[1:])
        
        fields = lines[0].strip().split("\t")  # Uses the fields in the database summary.
        if  self.fields ^ set(fields):  # Check whether the database file system and the code expect the same fields.
            if self.fields - set(fields):
                missing_fields = list(self.fields - set(fields))
                print("Warning: The field(s): %s not found in the database. Setting these attributes to None. Update the version of ABDB filesystem." % ", ".join(missing_fields), file=sys.stderr)
            else:
                print("Warning: Unimplemented field(s): %s found in the database. These will be ignored. Update the version of ABDB python package to use." % ", ".join(list(self.fields - set(fields))), file=sys.stderr)
        else:
            missing_fields = []
            
        for line in lines[1:]:
            values = line.strip().split("\t")
            if len(values) == len(fields):
                details = dict(list(zip(fields + missing_fields, values + [None] * len(missing_fields))))
                try:
                    self.db_summary[details["pdb"]]["fabs"].append(dict((fabf, details[fabf]) for fabf in self.fabfields))
                except KeyError:
                    self.db_summary[details["pdb"]] = dict((pdbf, details[pdbf]) for pdbf in self.pdbfields) 
                    self.db_summary[details["pdb"]]["fabs"] = [ dict((fabf, details[fabf]) for fabf in self.fabfields) ] 
                self.db_summary[details["pdb"]]["inhouse"]=inhouse

                    
    def _get_numbered_sequences(self, methods=["X-RAY DIFFRACTION"]):
        """
        Read the aligned sequences.
        Only load sequences for structures with method in methods.
        """
        if hasattr(self, "sequences"):
            if self.sequences:
                return self.sequences
        self.sequences = {}
        for p in self:
            pdb = self.fetch(p)
            if pdb.get_method() in methods:
                self.sequences[p] = {}
                for fab in pdb:
                    try:
                        self.sequences[p][fab] = fab.get_numbering()
                    except KeyError:
                        continue
        return self.sequences
    
    def _get_headerpath(self, pdb):
        """
        Get the path to the header file for the pdb.
        """
        headerpath = os.path.join(self.dbpath, "entries", pdb, "header", "%s.header.txt" % pdb)
        if os.path.exists(headerpath):
            return headerpath
        else:
            headerpath = os.path.join(self.dbpath, "inhouse_entries", pdb, "header", "%s.header.txt" % pdb)
            if os.path.exists(headerpath):
                return headerpath
    
    def _get_filepath(self, pdb):
        """
        Get the path to the structure file for the pdb.
        """
        filepath = os.path.join(self.dbpath, "entries", pdb, "structure", "%s.pdb" % pdb)
        if os.path.exists(filepath):
            return filepath
        else:
            filepath = os.path.join(self.dbpath, "inhouse_entries", pdb, "structure", "%s.pdb" % pdb)
            if os.path.exists(filepath):
                return filepath

    def _get_annotationpath(self, pdb):
        """
        Get the path to the annotation directory for the pdb.
        """
        annotationpath = os.path.join(self.dbpath, "entries", pdb, "annotation")
        if os.path.exists(annotationpath):
            return annotationpath
        else:
            annotationpath = os.path.join(self.dbpath, "inhouse_entries", pdb, "annotation")
            if os.path.exists(annotationpath):
                return annotationpath
    
    def _get_sequencepath(self, pdb):
        """
        Get the path to the sequence directory for the pdb.
        """
        sequencepath = os.path.join(self.dbpath, "entries", pdb, "sequences")
        if os.path.exists(sequencepath):
            return sequencepath
        else:
            annotationpath = os.path.join(self.dbpath, "inhouse_entries", pdb, "sequences")
            if os.path.exists(annotationpath):
                return annotationpath

    def _get_imgtpath(self, pdb):
        """
        Get the path to the imgt directory for the pdb.
        """
        imgtpath = os.path.join(self.dbpath, "entries", pdb, "imgt")
        if os.path.exists(imgtpath):
            return imgtpath
        else:
            imgtpath = os.path.join(self.dbpath, "inhouse_entries", pdb, "imgt")
            if os.path.exists(imgtpath):
                return imgtpath

    def _get_abanglepath(self, pdb):
        """
        Get the path to the abangle directory for the pdb.
        """
        abanglepath = os.path.join(self.dbpath, "entries", pdb, "abangle")
        if os.path.exists(abanglepath):
            return abanglepath
        else:
            abanglepath = os.path.join(self.dbpath, "inhouse_entries", pdb, "abangle")
            if os.path.exists(abanglepath):
                return abanglepath

    def _get_summarypath(self, pdb):
        """
        Get the path to the summary directory for the pdb.
        """
        summarypath = os.path.join(self.dbpath, "entries", pdb, "summary")
        if os.path.exists(summarypath):
            return summarypath
        else:
            summarypath = os.path.join(self.dbpath, "inhouse_entries", pdb, "summary")
            if os.path.exists(summarypath):
                return summarypath

    def _get_clusteringpath(self, pdb):
        """
        Get the path to the clustering directory
        """
        clusteringpath = os.path.join(self.dbpath, "clustering")
        if os.path.exists(clusteringpath):
            return clusteringpath
        
    def _validate_sequence(self, seq):
        """
        Check whether a sequence is a protein sequence or if someone has submitted something nasty.
        """
        if len(seq) > 1000:
            raise Exception("Sequence too long.")
        if any([1 for s in seq.upper() if s not in aa1]):
            raise Exception("Unknown amino acid letter found in sequence.")
        else:
            return True
        

class PDB_details(object):
    '''
    A class to hold the details about a pdb structure.
    This includes information about the antibody and antigen molecules in the structure.
    '''
    def __init__(self, id,database=None):
        self.id = id
        self.misc = {}
        self._database = database # Linking the PDB object to the database so that it can use database methods
        self._loaded_numbering_scheme = None 
        
    def __repr__(self):
        return self.id
    
    def __getitem__(self, x):
        """
        Get item can be used to return an Fab object from the object
        @param x: Heavy_ident+Light_ident
        @type x: str
        """
        for fab in self.get_fabs():
            if x == fab.VH + fab.VL:
                return fab
    
    def __iter__(self):
        for fab in self.get_fabs():
            yield fab
    
    # set values.
    def _set_method(self, method=None):
        self.method = method    

    def _set_resolution(self, resolution=None):
        self.resolution = resolution

    def _set_compound(self, compound=None):
        self.compound = compound.title()

    def _set_r_free(self, r_free=None):
        self.r_free = r_free

    def _set_r_factor(self, r_factor=None):
        self.r_factor = r_factor
    
    def _set_organism(self, organism=None):
        """
        Set the organism of the pdb - this is taken from the pdb.
        It is formatted in like:
        unique_organism1;  unique_organism2; unique_organism3
        if there are three different unique organisms in the pdb. 
        
        If a molecule is made of more than one organism (e.g. a chimera) the format would be:
        organism1a, organism1b;  unique_organism2; unique_organism3
        
        So we split the organism field into an internal variable and assign it to the antigen and fab molecules. (they can see _organisms)
        """
        self.organism = organism  # this is public and is the raw field.
        self._organisms = list(map(str.strip, organism.split(";")))  # this is private.
        
    def _set_species(self):
        """
        The species field is the species of the FABS in the pdb. This is NOT the same as the organism field which is the raw pdb field.
        In the rare case that there are fabs of different species in the same pdb, these two species will be separated by "; ".
        Chimeric fabs will be labelled "CHIMERIC species1/species2 /... "
        """ 
        self.species = "; ".join(uniq([ fab.get_species() for fab in self.get_fabs() ]))

    def _set_common_species(self):
        """
        The species field is the species of the FABS in the pdb. This is NOT the same as the organism field which is the raw pdb field.
        In the rare case that there are fabs of different species in the same pdb, these two species will be separated by "; ".
        Chimeric fabs will be labelled "CHIMERIC species1/species2 /... "
        """ 
        self.common_species = "; ".join(uniq([ fab.get_common_species() for fab in self.get_fabs() ]))

    def _set_antigen_species(self):
        self.antigen_species = "; ".join(uniq([ ag.get_species() for ag in self.get_antigens() ]))
         
    def _set_short_header(self, short_header=None):
        self.short_header = short_header
    
    def _set_authors(self, authors=None):
        self.authors = authors
    
    def _set_date(self, date=None):
        self.date = date

    def _set_filepath(self, filepath=None):
        """
        Set the filepath to the PDB structure.
        """
        self.filepath = filepath

    def _set_headerpath(self, headerpath=None):
        """
        Set the path to the header file for the structure
        """
        self.headerpath = headerpath
        
    def _set_annotationpath(self, annotationpath=None):
        """
        Set the path to the directory containing antibody numbering annotations.
        These annotations have been generated by abnum.
        http://www.bioinf.org.uk/abs/abnum/
        """
        self.annotationpath = annotationpath
        
    def _set_sequencepath(self, sequencepath=None):
        """
        Set the path to the directory containing sequences for the various chains in the structure.
        These are fasta files.
        """
        self.sequencepath = sequencepath

    def _set_imgtpath(self, imgtpath=None):
        """
        Set the path to the directory containing imgt gene annotations for the chains.
        """
        self.imgtpath = imgtpath

    def _set_abanglepath(self, abanglepath=None):
        """
        Set the path to the directory containing imgt gene annotations for the chains.
        """
        self.abanglepath = abanglepath

    def _set_summarypath(self, summarypath=None):
        """
        Set the path to the directory containing summary file for the entry
        """
        self.summarypath = summarypath      

    def _set_clusteringpath(self, clusteringpath=None):
        """
        Set the path to the directory containing information about the cdr clusters
        """
        self.clusteringpath = clusteringpath

    def _set_inhouse(self, inhouse=False ):
        """
        Set whether the structure is an inhouse structure or not.
        """            
        self.inhouse=inhouse
          
    def _set_fabs(self, fabs=None):
        """
        Set the fabs in the structure.
        fabs is a list of dictionaries containing information about each "fab".
        An fab might be a single chain if it has not been paired.
        """
        
        if fabs:
            self.fabs = []
            self.models = set()
            for fab in fabs:
                newfab = Fab_details(self, fab)
                self.fabs.append(newfab)
                self.models.add(newfab.model)
            if False in [fab.is_completefab() for fab in self.fabs]: 
                self._set_completefab(False)
            else:
                self._set_completefab(True)
            if any([fab.is_completefab() for fab in self.fabs]):
                self._set_contains_completefab(True)
            else:    
                self._set_contains_completefab(False)
                
            if any([fab.is_scfv() for fab in self.fabs]):
                self._set_scfv(True)
            else:
                self._set_scfv(False)
            if True in [fab.is_complex() for fab in self.fabs]:
                self._set_complex(True)
            else:
                self._set_complex(False)
                
            self.models = list(self.models)
        else:
            self.fabs = None
            self.models = []
            self._set_completefab(False)
            self._set_contains_completefab(False)
            self._set_scfv(False)
            self._set_complex(False)
            
    def _set_header(self):
        """
        Read the header from the database and set it as an attribute.
        """
        try:
            with open(self.headerpath) as header_file:
                self.header = header_file.read()
        except IOError:
            self.header = "Header not found on local database\nThis is the fake header for pdb %s\n" % self.id

    def _set_numbering(self, scheme=None ):
        """
        Set the antibody numbering for the chains as stored in the database.
        The numbering is done using the seqres sequence for each chain.

        
        Each annotation file describes the numbering of a variable region of a chain.
        Each chain may have more than one variable region (e.g. scfv)
        
        Annotations are stored for the pdb as:
        
        pdb_chain_region.ann
        
        e.g for 12e8 we have:
        12e8_H_VH.ann
        12e8_L_VL.ann
        12e8_M_VL.ann
        12e8_P_VH.ann
        
        These are Chothia numbered.
        
        self.numbering is a dictionary with chain identifiers as keys.
        
        For each chain there is a list of variable regions.
        
        Each entry of the list is a tuple of chain_numbering and chain_type.
        
        chain_type will be 'H' or 'L'.
        
        chain numbering is a list of biopython.PDB residue identifiers and one letter amino acid codes.
        
        e.g. [ ( (1,' '), E ), ( (2, ' '), Q ) ..... ]
        
        This is not pretty but it works... 
        We may want to re-think this in the future.
        
        """
        def interpret(x):
            """
            Conversion from abnum residue identifiers to biopython.PDB residue identifiers.
            e.g 'H100A'  --> ( 100, 'A') 
            e.g 'H100'   --> ( 100, ' ')
            """
            assert x[0] == "H" or x[0] == "L", x
            try:
                return (int(x[1:]), ' ')
            except ValueError:
                return (int(x[1:-1]), x[-1])
        
        try:
            if self.annotationpath is None:
                raise OSError

            if scheme is None:
                scheme=self._database.get_numbering_scheme()

            if scheme == "chothia": # Back compatibility - to be removed in time
                if os.path.exists( os.path.join( self.annotationpath, scheme ) ): # If new system use new file location
                    usepath = os.path.join( self.annotationpath, scheme )
                else:
                    usepath = self.annotationpath # If old system use old file location
            else: 
                assert os.path.exists( os.path.join( self.annotationpath, scheme ) ), "%s scheme does not have a numbering file for %s"%(scheme, self.id) # Ensure the directory exists

            usepath = os.path.join( self.annotationpath, scheme )

            self.numbering = {}

            annotations = os.listdir(usepath)

            for fname in annotations:
                root, ext = os.path.splitext(fname)
                if ext == ".ann":
                    pdb_name, chain, region = root.split("_")
                    if pdb_name == self.id:
                        with open(os.path.join(usepath, fname)) as f:
                            annotation = f.read() 
                            if not annotation:
                                self._loaded_numbering_scheme = None
                                return

                            result = list(map(str.split, annotation.strip().split("\n")))
                            chain_numbering = [ (interpret(res[0]), res[1])  for res in result if res[1] != "-" ]
                            chain_type = result[0][0][0]
                            assert chain_type == region[1]
                        try:
                            self.numbering[chain].append((chain_numbering, chain_type))
                        except  KeyError:
                            self.numbering[chain] = [(chain_numbering, chain_type)]
            self._loaded_numbering_scheme = scheme
        except OSError:
            self._loaded_numbering_scheme = None
            pass
    
    def _set_alignment(self, scheme=None):
        """
        An alignment is needed between the seqres sequence and the structural sequence (missing residues in the structure will lead to differences).

        This retrieves the alignment from the database
        
        We create an alignment dictionary between the seqres sequence and the structure sequence. 
        
        See ABDB.Annotate.annotate.get_alignment_dict for more details.
        """

        try:
            # Check we can reach the sequences
            if self.sequencepath is None:
                raise OSError

            if scheme is None:
                scheme = self._database.get_numbering_scheme() # Always default to the known. 

            if scheme == "chothia": # Back compatibility - to be removed in time
                if os.path.exists( os.path.join( self.sequencepath, scheme ) ): # If new system use new file location
                    usepath = os.path.join( self.sequencepath, scheme )
                else:
                    usepath = self.sequencepath # If old system use old file location
            else: 
                assert os.path.exists( os.path.join( self.sequencepath, scheme ) ), "%s scheme does not have an alignment file for %s"%(scheme, self.id) # Ensure the directory exists

            usepath = os.path.join( self.sequencepath, scheme )

            self.alignment = {} 

            numbering = self.get_numbering(scheme=scheme) # Gets the numbering required

            self.sequence = {}
            for chain in numbering:
                self.sequence[chain] = {}
                for region in numbering[chain]:
                    with open(os.path.join(usepath, "%s_%s_V%s.fa" % (self.id, chain, region[1]))) as f:
                        s = f.read().split(">")
                        self.sequence[chain][region[-1]] = {}
                        for e in s:
                            if e:
                                lines = e.split("\n")
                                self.sequence[chain][region[-1]][ "".join(lines[0].split("|")[1:]) ] = lines[1]
                        # Get the alignment between the structure and numbering 
                        try:
                            self.alignment[chain].append(get_alignment_dict(self.sequence[chain][region[-1]]["structurefull"], self.sequence[chain][region[-1]]["seqresregion: %s" % region[1]]))
                        except KeyError:
                            self.alignment[chain] = [ get_alignment_dict(self.sequence[chain][region[-1]]["structurefull"], self.sequence[chain][region[-1]]["seqresregion: %s" % region[1]]) ]

        except OSError:
            pass

    def _set_sequence(self, scheme=None):
        self._set_alignment(scheme=scheme)

    def _set_raw_sequences(self): # All the sequences in the file
        self.raw_sequences = {}
        try:
            with open(os.path.join(self.sequencepath,"%s_raw.fa"%self.id) ,'r' ) as f:
                for line in f:
                    stripped = line.strip()
                    if stripped[0]=="#":
                        continue
                    elif stripped[0] == ">":
                        chain = stripped.split("|")[0][-1]
                        self.raw_sequences[chain]=""
                    elif line.strip():
                        self.raw_sequences[chain] += stripped
        except (IOError, KeyError) as e:
            self.raw_sequences = {}
                
      
    def _set_antigens(self, antigens=None):
        if antigens:
            self.antigens = []
            for antigen in antigens:
                self.antigens.append(antigen)
        else:
            self.antigens = []
        
    def _set_complex(self, complex=None):
        self.complex = complex
    
    def _set_scfv(self, scfv=None):
        self.scfv = scfv      
    
    def _set_completefab(self, completefab=None):
        self.completefab = completefab

    def _set_contains_completefab(self, contains_completefab=None):
        self.contains_completefab = contains_completefab

    def _set_engineered(self, engineered=None):
        self.engineered = engineered      
      
    def _set_affinity(self, affinity=None):
        """
        Affinity data from Jin's database.
        
        @param affinity: The affinity data collected (Kd) - or None
        @type affinity: C{float} or C{str} or None. Expect an exponent format for the string or float. e.g 1.0e-10 
        """
        try:
            self.affinity = float(affinity)
        except ValueError:  # Got "None" /NA etc
            self.affinity = None
        except TypeError:  # Got None
            self.affinity = None

    def _set_delta_g(self, delta_g=None):
        """
        Delta G data from Jin's database.
        
        @param delt_g: The calculated delta G value
        @type delta g: C{float} or C{str} or None.  
        """
        try:
            self.delta_g = float(delta_g)
        except ValueError:  # Got "None" /NA etc
            self.delta_g = None
        except TypeError:  # Got None
            self.delta_g = None

    def _set_pmid(self, pmid=None):
        """
        Pubmed id of the relevent publication - currently only for structures with affinity data
        
        @param pmid: The pubmed id of the paper
        """
        try:
            int(pmid)
            self.pmid = pmid
        except ValueError:  # Got "None" /NA etc
            self.pmid = None
        except TypeError:  # Got None
            self.pmid = None
            
    def _set_affinity_method(self, affinity_method=None):
        """
        """
        if affinity_method is not None:
            if affinity_method.lower() == "none" or affinity_method.lower() == "na":
               affinity_method=None
        self.affinity_method=affinity_method

    def _set_temperature(self, temperature=None):
        """
        """
        try:
            float(temperature)
            self.temperature = temperature
        except ValueError:  # Got "None" /NA etc
            self.temperature = None
        except TypeError:  # Got None
            self.temperature = None

    
    def _set_imgt_details(self):
        """
        Set the imgt details of all the antibody chains in the pdb.
        """
        
        self.imgt_details = {}
        
        try:
            if self.imgtpath is None:
                raise OSError
            imgt_files = os.listdir(self.imgtpath)
            for fname in imgt_files:
                try:
                    pdb, chain_id, chain_type = os.path.splitext(fname)[0].split("_")
                    
                    with open(os.path.join(self.imgtpath, fname)) as f:
                        lines = f.readlines()
                        header = lines[0].split()
                        self.imgt_details[ (pdb, chain_id, chain_type) ] = {}
                        for line in lines[1:]:
                            if line.strip():
                                details = dict(list(zip(header, line.split())))
                                self.imgt_details[ (pdb, chain_id, chain_type) ][details["gene_type"]] = details
                except ValueError:
                    # fname was incorrect format.
                    continue
                except IOError:
                    # could not open the file.
                    continue
                except IndexError:
                    # The file was empty
                    continue
                except KeyError:
                    # The header did not have the field "gene_type"
                    continue
        except OSError:
            # The imgtpath was not set
            pass

    def _set_orientation_angles(self):
        """
        Read the orientation angles from the database.
        """
        self.orientation_angles = {}
        if self.abanglepath:
            try:
                with open(os.path.join(self.abanglepath, self.id + ".abangle")) as f:
                    lines = f.readlines()
                    if lines:
                        header = lines[0].strip().split("\t")
                        for line in lines[1:]:
                            details = dict(list(zip(header, line.strip().split("\t"))))
                            self.orientation_angles[(details["pdb"], details["Hchain"], details["Lchain"], int(details["model"]))] = dict((a, float(details[a])) for a in ["HL", "HC1", "HC2", "LC1", "LC2", "dc"])
                            
            except IOError:
                pass
            except IndexError:
                pass
                
    def _set_summary_lines(self):
        """
        Read the summary lines from the database.
        """
        self.summary_lines = {}
        if self.summarypath:
            try:
                with open(os.path.join(self.summarypath, self.id + ".tsv")) as f:
                    lines = f.readlines()
            except IOError:
                try:
                    with open(os.path.join(self.summarypath, self.id + ".csv")) as f:
                        lines = f.readlines()
                except IOError:
                    lines = []
            if lines:
                self.summary_lines["header"] = lines[0].strip()
                for line in lines[1:]:
                    pdb,H,L,model=line.split("\t")[:4]
                    if H.upper()==L.upper(): # deal with structural scfvs
                        H=H.upper()
                        L=L.upper()
                        self.summary_lines[tuple([pdb,H,L,model])]="\t".join(line.split("\t")[:1] + list(map(str.upper, line.split("\t")[1:3])) + line.strip().split("\t")[3:])  
                    else:
                        self.summary_lines[tuple([pdb,H,L,model])]=line.strip()
    
    def _load_structure(self,scheme=None,definition=None):
        """
        Load the structure from the Database.
        This sets the structure attribute.
        Use get_structure to return the structure. 
        """
        try:
            with open( os.path.join(os.path.split(self.filepath)[0], "ignore.dat" ) ) as ignorefile:
                ignore_contacts = [tuple(l.strip().split()) for l in ignorefile.readlines()]
        except IOError:
            ignore_contacts=[]
          
        # Get the numbering scheme to be used.
        if scheme is None: # If it is None then use the the default scheme 
            scheme = self._database.get_numbering_scheme()
        if definition is None:
            definition = self._database.get_region_definition()

        numbering = self.get_numbering(scheme)
        parser = AntibodyParser(QUIET=True,scheme=scheme,definition=definition)
        self.structure = parser.get_antibody_structure(self.id, self.get_filepath(), numbering, self.get_alignment(),crystal_contacts=ignore_contacts)
      
    def clear_structure(self):
        """
        This deletes the structure attribute if it has been loaded.
        Note: if you have other references to the structure object they MUST ALSO BE DELETED to free up the memory.
        e.g. s = pdb_details.get_structure()
             pdb_details.clear_structure()
             del s
             
        For batch analysis of structures with this code this should be considered.
        No memory considerations are made (i.e. do not load in ALL the structures at the same time or your machine will die :( )
        """
        if hasattr(self, "structure"):
            del self.structure

    # Get Methods
    def get_method(self):
        """
        Get the experimental method used to obtain the structure.
        """
        return self.method    

    def get_resolution(self):
        """
        Get the resolution of the structure. 
        This only applies to X-ray structures.
        """
        return self.resolution
    
    def get_r_factor(self):
        """
        Get the r_factor (observed) of the structure. 
        This only applies to X-ray structures.
        """
        return self.r_factor
    
    def get_r_free(self):
        """
        Get the r_free of the structure. 
        This only applies to X-ray structures.
        """
        return self.r_free
    
    def get_header(self):
        """
        Get the header of the PDB file
        """
        if not hasattr(self, "header"):
            self._set_header()
        elif not self.header:
            self._set_header()

        return self.header

    def get_short_header(self):
        """
        Get the short header of the PDB file.
        This often contains the molecule name or title.
        """
        return self.short_header
    
    def get_authors(self):
        """
        Get the authors of the structure.
        """
        return self.authors

    def get_date(self):
        """
        Get the date deposited in the PDB.
        """
        return self.date

    def get_pmid(self):
        """
        Get the pubmed id of the structure.
        
        Currently only returns pmid for structures with affinity data.
        """
        return self.pmid

    def get_compound(self):
        """
        Get the name of the compound
        """
        return self.compound

    def get_organism(self):
        """
        Get the organism(s) the molecule has been obtained from.
        Note that this is the raw organism field from the PDB. NOT the same as get_species.
        Therefore, those structures that are more than one species will have both species names in the field.
        """
        return self.organism
    
    def get_species(self):
        """
        Get the species of the antibodies chains in the molecule. 
        If the fabs have different species (rare) then the seperate species will be returned separated with a ";".
        Chimeric antibodies should have two species separated by a ",".
        """
        if not hasattr(self, "species"):
            self._set_species()
        return self.species    

    def get_common_species(self):
        """
        Get the species of the antibodies chains in the molecule. 
        If the fabs have different species (rare) then the seperate species will be returned separated with a ";".
        Chimeric antibodies should have two species separated by a ",".
        """
        if not hasattr(self, "common_species"):
            self._set_common_species()
        return self.common_species    
    
    def get_antigen_species(self):
        """
        Get the species of the antigen molecule if present. 
        This is only defined for biological protein and peptide antigens. 
        """
        if not hasattr(self, "antigen_species"):
            self._set_antigen_species()
        return self.antigen_species    
    
    def get_filepath(self):
        """
        Get the path to the PDB structure in the database.
        """
        return self.filepath

    def get_fabs(self):
        """
        Get a list of the fabs (or unpaired antibody chains) in the PDB structure.
        Each item of the list will be an fab_details object.
        """
        return self.fabs

    def get_numbering(self, scheme=None ):
        """
        Get the numbering of the chains in the PDB file.
        See L{set_numbering<self._set_numbering>} for more information about the numbering attribute.
        """
        if scheme is None: # Then get the default scheme for the database
            scheme = self._database.get_numbering_scheme()
        if self._loaded_numbering_scheme != scheme:
            self._set_numbering(scheme=scheme)
        elif not hasattr(self, "numbering"):
            self._set_numbering(scheme=scheme)
        elif not self.numbering:
            self._set_numbering(scheme=scheme)

        return self.numbering
    
    def get_alignment(self, scheme=None):
        """
        Get the alignment dictionary between the seqres sequence and the structure sequence.

        This is a dictionary with chain IDs for keys and alignment dictionaries for values.
        
        Alignment dictionary tells you that the the i'th (key) position in the seqres sequence corresponds to the j'th (value) position in the structure sequence.
        
        See ABDB.Annotate.annotate.get_alignment_dict for more details.
        """

        if scheme is None: # Then get the default scheme for the database
            scheme = self._database.get_numbering_scheme()
        if self._loaded_numbering_scheme != scheme: # If the user has requested a different scheme
            self._set_alignment(scheme=scheme)
        elif not hasattr(self, "alignment"):
            self._set_alignment(scheme=scheme)
        elif not self.numbering:
            self._set_alignment(scheme=scheme)

        return self.alignment
    
    def get_structure(self, scheme=None,definition=None):
        """
        Get the structure loaded as a biopython pdb structure put into an "antibody context".
        The antibody chains will have been renumbered using the Chothia numbering system. See the ABDB.AbPDB package for more information.
        """

        if scheme is None: # Then get the default scheme for the database
            scheme = self._database.get_numbering_scheme()
        if definition is None:
            definition = self._database.get_region_definition()

        if self._loaded_numbering_scheme != scheme: # If the user has requested a different scheme
            self._load_structure(scheme=scheme, definition=definition)
        elif not hasattr(self, "structure"):
            self._load_structure(scheme=scheme, definition=definition)
        elif not self.numbering:
            self._load_structure(scheme=scheme, definition=definition)

        return self.structure


    def get_sequence(self, scheme=None):
        """
        Get the sequences of the antibody chains in the structure.
        
        Returns a nested dictionary.
        
        For example a PDB has antibody chains P and M. P is a heavy chain, M is a light chain.
        
        s = pdb_object.get_sequence()
        
        s.keys() # [ "P" , "M" ]
        
        s["P"] # This returns a dictionary for each variable region in the chain ( single chain fvs can have more than one variable region e.g. 2kh2).         
        
        s["P"]["H"] # The heavy region of chain P          
        s["M"]["L"] # The light region of chain M
        
        s["P"]["H"] # contains:
                    #    structurefull    - the sequence found in the structure (i.e. resolved residues)
                    #    seqresfull       - the sequence found in the seqres record from the pdb
                    #    seqresregion: H  - the variable region of the seqres record ( non-variable region is represented by "-" )
                    
        s["M"]["L"] # contains:
                    #    structurefull    - the sequence found in the structure (i.e. resolved residues)
                    #    seqresfull       - the sequence found in the seqres record from the pdb
                    #    seqresregion: L  - the variable region of the seqres record ( non-variable region is represented by "-" )
          
        """


        if scheme is None: # Then get the default scheme for the database
            scheme = self._database.get_numbering_scheme()
        if self._loaded_numbering_scheme != scheme: # If the user has requested a different scheme
            self._set_sequence(scheme=scheme)
        elif not hasattr(self, "sequence"):
            self._set_sequence(scheme=scheme)
        elif not self.sequence:
            self._set_sequence(scheme=scheme)

        return self.sequence

    def get_raw_sequences(self):
        """
        Get a dictionary containing the chain identifier as the key and the sequence as the value.
        """
        if not hasattr(self, "raw_sequence"):
            self._set_raw_sequences()
        elif not self.raw_sequences:
            self._set_raw_sequences()

        return self.raw_sequences

    def get_summary_lines(self):
        """
        Get the summary lines for the structure.
        This will return a dictionary with either 'header' or a fab identifier ( i.e (pdb, H, L, modelid) ) as keys and the line as a value. 
        """
        if not hasattr(self, "summary_lines"):
            self._set_summary_lines()
        elif not self.summary_lines:
            self._set_summary_lines()
        return self.summary_lines

    def get_antigens(self):
        """
        Get a list of the antigen molecules in the PDB file.
        Each item of the list will be an antigen_details object.
        If there are no antigens this will be an empty list.
        """
        return self.antigens
    
    def get_light_chain_types(self):
        """
        Get a list of the light chain types of the fabs in the pdb (normally these will all be the same)
        """
        return [ fab.get_light_chain_type() for fab in self.get_fabs()]
    
    def get_light_subclasses(self):
        """
        Get a list of the light subclass of the fabs in the pdb (normally these will all be the same)
        """
        return [ fab.get_light_subclass() for fab in self.get_fabs()]
    
    def get_heavy_subclasses(self):
        """
        Get a list of the light subclass of the fabs in the pdb (normally these will all be the same)
        """
        return [ fab.get_heavy_subclass() for fab in self.get_fabs()]
    
    def get_imgt_details(self):
        """
        Get a dictionary containing the imgt details for all the antibody chains in this pdb.
        Query individual fabs for their details.
        """
        if not hasattr(self, "imgt_details"):
            self._set_imgt_details()
        elif not self.imgt_details:
            self._set_imgt_details()

        return self.imgt_details
    
    def get_orientation_angles(self):
        """
        Get a dictionary containing the abangle angles for each fv in this pdb.
        """
        if not hasattr(self, "orientation_angles"):
            self._set_orientation_angles()
        elif not self.orientation_angles:
            self._set_orientation_angles()

        return self.orientation_angles

    def is_complex(self):
        """
        Check whether the Fabs in the structure are bound to antigens or not. 
        This is true if any of the fabs are bound and false otherwise.
        Inspect each fab individually to check their bound state
        """
        return self.complex
    
    
    def get_affinity(self):
        """
        Get the affinity (Kd) of the antibody in this structure to its antigen.
                
        @return: The affinity (Kd) of the antibody to its antigen. Units are M. 
        @rtype: C{float} or None
        """
        return self.affinity

    def get_delta_g(self):
        """
        Get the delta G value for the antibody in this structure to its antigen.
                
        @return: The delta G value of the antibody to its antigen. Units are kcal/mol 
        @rtype: C{float} or None
        """
        return self.delta_g
    
    def get_temperature(self):
        """
        Get the temperature at which the affinity information was obtained.
        """
        return self.temperature

    def get_affinity_method(self):
        """
        Get the method by which the affinity information was obtained.
        """
        return self.affinity_method

    
    def is_scfv(self):
        """
        Check whether any of the Fabs in the structure are single chain fvs.
        This is parsed from the header record of the PDB using nlp. 
        If their are multiple variable regions on the chain these will also be a scfv.
        """
        return self.scfv
    
    def is_completefab(self):
        """
        Check whether the fvs in the structure are complete i.e. paired. 
        If there are any unpaired antibody chains within the structure this returns False, otherwise True.
        * set to be depreciated as is misleading - in favour of has_completefab method
        """
        return self.completefab

    def has_completefab(self):
        """
        Check whether the pdb contains a paired VH-VL. 
        Returns True if it does False if not.
        Inspect each individual fab details object to see which ones are paired.
        
        Slightly misleading. 
        """
        return self.contains_completefab

        
    def is_engineered(self):
        """
        Check whether the structure is engineered or not.
        This has been parsed from the header.
        """
        return self.engineered
    
    def is_imgt_from_imgt(self):
        """
        Check whether source of the imgt details is from the imgt webserver.
        """
        imgt_details = self.get_imgt_details()
        if imgt_details:
            if list(imgt_details.values())[0].values()[0]["parsed_from_imgt"] == "True":
                return True
            else:
                return False
        else:
            return False

    def is_imgt_from_ABDB(self):
        """
        Check whether source of the imgt details is assigned by ABDB.
        """
        imgt_details = self.get_imgt_details()
        if imgt_details:
            if list(imgt_details.values())[0].values()[0]["parsed_from_imgt"] == "False":
                return True
            else:
                return False
        else:
            return False

    def is_inhouse(self):
        """
        Is the entry an in-house structure?
        """        
        return self.inhouse


class Fab_details(object):
    '''
    Class to hold the details of an fab
    '''

    def __init__(self, pdb, fab):
        self._set_pdb(pdb)
        self._set_chains(fab["Hchain"] , fab["Lchain"])        
        self._set_model(fab["model"])
        self.id = (self.pdb.id, self.VH, self.VL, int(self.model))  # for consistency with structure fab id
        self._set_heavy_subclass(fab["heavy_subclass"])  # these should be change to "variable" subclass. Further gene annotation can then be read in from imgt
        self._set_light_subclass(fab["light_subclass"])
        self._set_light_chain_type(fab["light_ctype"])
        self._set_species(fab["heavy_species"], fab["light_species"])
        self._set_scfv(fab["scfv"])
        self._set_completefab()
        
        if fab["antigen_chain"] == "NA":
            antigen = None
        else:
            antigen = Antigen_details(fab["antigen_chain"]) 
            antigen._set_pdb(pdb)
            antigen._set_antigen_type(fab["antigen_type"]) 
            antigen._set_antigen_het_name(fab["antigen_het_name"])  # this is added functionality - NA if polymer. hetname if small molecule
            antigen._set_antigen_name(fab["antigen_name"])  # this is added functionality - name of the antigen parsed from chemical component db or header.
            antigen._set_species(fab["antigen_species"])  # this is added functionality - name of the antigen parsed from chemical component db or header.
            antigen._set_fab(self)

        self._set_antigen(antigen)
        self.numbering = None
        self.sequence = None

    def __str__(self):
        if self.completefab:
            return "fab %s%s\n\tVH = chain %s\n\tVL = chain %s" % (self.VH, self.VL, self.VH, self.VL)
        elif self.VH == "NA":
            return "ab_chain %s\n\tVL = chain %s\n" % (self.VL, self.VL)
        elif self.VL == "NA":
            return "ab_chain %s\n\tVH = chain %s\n" % (self.VH, self.VH)

    def __repr__(self):
        if self.completefab:
            return "\nfab %s%s\n\tVH = chain %s\n\tVL = chain %s\n" % (self.VH, self.VL, self.VH, self.VL)
        elif self.VH == "NA":
            return "\nab_chain %s\n\tVL = chain %s\n" % (self.VL, self.VL)
        elif self.VL == "NA":
            return "\nab_chain %s\n\tVH = chain %s\n" % (self.VH, self.VH)
        
    def _set_pdb(self, pdb):
        self.pdb = pdb  # *
        
    def _set_chains(self, h, l):
        """
        Deal with being a structural scfv - it will come in as Bb or bB instead of BB.
        """
        if l.upper() == h.upper():        
            if l in self.pdb.get_sequence():  # the lower case is actually correct
                self._set_VL(l)
            else:
                self._set_VL(l.upper())
            if h in self.pdb.get_sequence():  # the lower case is actually correct
                self._set_VH(h)
            else:
                self._set_VH(h.upper())
        else:
            self._set_VL(l)
            self._set_VH(h)
    
    def _set_VH(self, vh):
        self.VH = vh

    def _set_VL(self, vl):
        self.VL = vl

    def _set_model(self, model):
        self.model = model
    
    def _set_heavy_subclass(self, subclass=""):
        self.heavy_subclass = subclass

    def _set_light_subclass(self, subclass=""):
        self.light_subclass = subclass
    
    def _set_light_chain_type(self, chaintype=""):
        self.light_chain_type = chaintype
        
    def _set_species(self, hspecies=None, lspecies=None):
        """
        Set the species of the fab. 
        """
        if hspecies and lspecies:  # we have organisms for both 
            if hspecies == lspecies:  # they are the same - easier.
                    if "," in hspecies:
                        self.species = "CHIMERIC %s" % "/".join(sorted(map(str.strip, hspecies.split(",")))).upper()
                        # I sort the above so that it always gives the same ordering of species - they are not in any standardized order anyway.
                    else:
                        self.species = hspecies.upper()
            else:
                s = hspecies + "," + lspecies
                self.species = "CHIMERIC %s" % "/".join(sorted(uniq(list(map(str.strip, s.split(",")))))).upper()
        elif hspecies:
            if "," in hspecies:
                self.species = "CHIMERIC %s" % "/".join(sorted(map(str.strip, hspecies.split(",")))).upper()
                # I sort the above so that it always gives the same ordering of species - they are not in any standardized order anyway.
            else:
                self.species = hspecies.upper()
        elif lspecies:
            if "," in lspecies:
                self.species = "CHIMERIC %s" % "/".join(sorted(map(str.strip, lspecies.split(",")))).upper()
                # I sort the above so that it always gives the same ordering of species - they are not in any standardized order anyway.
            else:
                self.species = lspecies.upper()
                # 1. Check the length of pdb._organisms. If there is only one entry make the species equal to that entry
        elif len(self.pdb._organisms) == 1:
            # check for chimeric.
            if "," in self.pdb._organisms[0]:
                self.species = "CHIMERIC %s" % "/".join(sorted(map(str.strip, self.pdb._organisms[0].split(",")))).upper()
                # I sort the above so that it always gives the same ordering of species - they are not in any standardized order anyway.
            else:
                self.species = self.pdb._organisms[0].upper()
        else:
            self.species = "SYNTHETIC CONSTRUCT"
        
        untrusted_ab_species = ["ESCHERICHIA COLI", "CLOSTRIDIUM BOTULINUM"]  # add as necessery.
        if self.species == "LAMA": self.species = "LAMA GLAMA"  # default Lama
        if self.species == "SYNTHETIC": self.species = "SYNTHETIC CONSTRUCT"  # default synthetic
        if self.species == "UNCLASSIFIED": self.species = "UNKNOWN"  # default unknown
        if self.species in  untrusted_ab_species: self.species = "UNKNOWN"       
        if not self.species: self.species = "UNKNOWN"
        self._set_common_species()
        
    def _set_common_species(self):
        try:
            self.common_species = sci_to_common_names[self.get_species().lower()].upper()
        except KeyError:
            self.common_species = self.species
        
    def _set_complex(self, complex=False):
        """
        Boolean in complex with an antigen
        """
        self.complex = complex
    
    def _set_antigen(self, antigen=None):
        """
        Link to antigen
        """
        if antigen:
            self._set_complex(True)
            self.antigen = antigen
        else:
            self._set_complex(False) 
            self.antigen = None
    
    def _set_scfv(self, scfv=None):
        if scfv == "True":
            self.scfv = True
        elif scfv == "False":
            self.scfv = False

    def _set_completefab(self):
        if self.VH == "NA" or self.VL == "NA":
            self.completefab = False
        else:
            self.completefab = True

    def _set_CDR_sequences(self,scheme=None):
        """
        Set the cdr sequences using the Chothia, Kabat, Contact, IMGT and North defintions
        The residues in the seqres sequence are used (not the structure which may have missing residues)        

        @change: Now using the new Accept object from regions module and compatibility with the different schemes
        @change: Added heuristics to do conversion of non-equivalent positions.
        @change: Now use the annotate regions function as a utility function from AB_Utils
        """

        chains  = list(self.get_numbering(scheme=scheme).keys())

        cids = { "H": self.get_VH(), "L": self.get_VL() }

        self.CDR_sequences = {"chothia":{}, "kabat":{}, "contact":{}, 'imgt':{}, "north":{}}
        for chain in chains:
            numbering = [ _[0]  for _ in self.pdb.get_numbering(scheme=scheme)[ cids[chain] ] if _[1] == chain ][0]
            for definition in self.CDR_sequences:
                regions = annotate_regions( numbering, chain, numbering_scheme=scheme, definition=definition )
                for cdr in [1,2,3]:
                    cdr_str = ("cdr%s%d"%(chain,cdr)).lower()
                    self.CDR_sequences[definition][ cdr_str.upper() ] = [ (res,a) for res,a,r in regions if r==cdr_str ]

    def _set_CDR_lengths(self,scheme=None):
        """
        Set the cdr lengths using the Chothia, Kabat and Contact definitions.
        The residues in the seqres sequence are used (not the structure which may have missing residues)        
        """
        sequences = self.get_CDR_sequences(definition="*", scheme=scheme)
        self.CDR_lengths = {"chothia":{}, "kabat":{}, "contact":{}, "imgt":{}, "north":{}}
        for definition in ["chothia", "kabat", "contact", "imgt", "north"]:
            for cdr in sequences[definition]:
                self.CDR_lengths[definition][cdr] = len(sequences[definition][cdr])
                
                
    def _set_orientation_angles(self):
        """
        Set the orientation angles between the heavy and light variable domain for this fab.
        These are the ABangle angles.
        """
        try:
            self.orientation_angles = self.pdb.get_orientation_angles()[self.id]
        except KeyError:
            if self.VH == self.VL: # Fix for SCFV's to call orientation_angles from PDB_Details
                self.orientation_angles = dict([ ("".join(k[1:-1]).upper(), v) for (k,v) in list(self.pdb.get_orientation_angles().items()) ])[ "".join(self.id[1:-1]) ]
            else:
                self.orientation_angles = {}

    def _set_numbering(self,scheme=None):
        """
        Get the numbering from the parent pdb details object.
        """

        pdb_numbering = self.pdb.get_numbering(scheme=scheme)
        if not pdb_numbering:
            if self.VL == "NA":
                self.numbering = {"H": {}}
            elif self.VH == "NA":
                self.numbering = {"L": {}}
            else:
                self.numbering = {"H": {}, "L": {}}

        else:
            if self.VH==self.VL and len(pdb_numbering[self.VH])==2:
                self.numbering=dict( (x[1], dict(x[0])) for x in pdb_numbering[self.VH] ) # gets the H and L
            elif self.completefab:
                self.numbering = {"H":dict(pdb_numbering[self.VH][0][0]), "L":dict(pdb_numbering[self.VL][0][0])}
            elif self.VH == "NA":
                self.numbering = {"L":dict(pdb_numbering[self.VL][0][0])}
            elif self.VL == "NA":
                self.numbering = {"H":dict(pdb_numbering[self.VH][0][0])}
        
    def _set_sequence(self,scheme=None):
        """
        Set the sequence of the fv from the parent pdb.
        """
        pdb_sequence =    self.pdb.get_sequence(scheme=scheme)

        try:
            if self.VH == "NA":
                self.sequence = "-" + "/" + pdb_sequence[self.VL]["L"]["seqresregion: L"].strip("-")
            elif self.VL == "NA":
                self.sequence = pdb_sequence[self.VH]["H"]["seqresregion: H"].strip("-") + "/" + "-"
            else:
                self.sequence = pdb_sequence[self.VH]["H"]["seqresregion: H"].strip("-") + "/" + pdb_sequence[self.VL]["L"]["seqresregion: L"].strip("-")
        except KeyError: # Handles database corruption if not sequence file is found. 
            self.sequence = ""

        
    def _set_imgt_details(self):
        """
        Set the imgt details for the fab.
        """    
        pdb_imgt_details = self.pdb.get_imgt_details()
        
        if self.VH == "NA":
            self.heavy_imgt_details = {}
        else:
            try:
                self.heavy_imgt_details = pdb_imgt_details[(self.pdb.id, self.VH, "H")]
            except KeyError:
                self.heavy_imgt_details = {}
        
        if self.VL == "NA":
            self.light_imgt_details = {}
        else:
            try:
                self.light_imgt_details = pdb_imgt_details[(self.pdb.id, self.VL, "L")]
            except KeyError:
                self.light_imgt_details = {}
        
    def _set_missing_residues(self, scheme=None):
        self.missing_residues={}    
        if self.VH!="NA":
            seqres=self.pdb.get_sequence(scheme=scheme)[self.VH]["H"]["seqresregion: H"].rstrip()
            struc=self.pdb.get_sequence(scheme=scheme)[self.VH]["H"]["structurefull"]
            i,missing=0,[]
            for rseq, rstruc in zip(seqres, struc):
                if rstruc=="-"and rseq!="-": missing.append(i)
                if rseq!="-":i+=1
            if missing:
                if self.VH==self.VL:
                    i=[x[1] for x in self.pdb.get_numbering(scheme=scheme)[self.VH]].index("H")
                else:
                    i=0
                self.missing_residues["H"]=[ self.pdb.get_numbering(scheme=scheme)[self.VH][i][0][m][0] for m in missing ]
        if self.VL!="NA":
            seqres=self.pdb.get_sequence(scheme=scheme)[self.VL]["L"]["seqresregion: L"].rstrip()
            struc=self.pdb.get_sequence(scheme=scheme)[self.VL]["L"]["structurefull"]
            i,missing=0,[]
            for rseq, rstruc in zip(seqres, struc):
                if rstruc=="-"and rseq!="-": missing.append(i)
                if rseq!="-":i+=1
            if missing:
                if self.VH==self.VL:
                    i=[x[1] for x in self.pdb.get_numbering(scheme=scheme)[self.VL]].index("L")
                else:
                    i=0
                self.missing_residues["L"]=[ self.pdb.get_numbering(scheme=scheme)[self.VL][i][0][m][0] for m in missing ]

    # JD 13-03-15
    def _set_identicals(self, threshold=0.99):
        """
        Get all fab id's from the database that are similar by sequence id.
        """
        # No numbering scheme allowed - must use the database default as it will call the _get_sequences function.
        get_fabs = False
        if self.VH==self.VL or self.completefab:
            get_fabs = True
        # Filter out this fab
        identicals_and_me = self.get_pdb()._database.get_identical(self.get_numbering(), threshold=threshold, get_fabs=get_fabs)
        self.identicals = dict( _ for _ in identicals_and_me.items() if self.id != _[0] )
        self._loaded_threshold = "%.2f"%threshold # string comparison for reload check
        
    # JL 12-03-15
    # JD 13-03-15 mod
    def _set_unbound_forms(self,threshold=0.99):
        """
        Set the other unbound forms of this fv from the database
        """
        self.unbound_forms = []
        # There might be multiple unbounds.
        identicals=self.get_identicals( threshold=threshold )
        for iden, fab in identicals.items():
            if fab.antigen:
                continue
            else:
                self.unbound_forms.append( fab )
                
    # JL 12-03-15
    # JD 13-03-15 mod
    def _set_bound_forms(self,threshold=0.99):
        """
        Set the other bound forms of this fv from the database
        """
        self.bound_forms = []
        # There might be multiple bounds.
        identicals=self.get_identicals( threshold=threshold )
        for iden, fab in identicals.items():
            if not fab.antigen:
                continue
            else:
                self.bound_forms.append( fab )

    def get_pdb(self):
        """
        Get the pdb_details object that the fab belongs to.
        """
        return self.pdb
    
    def get_parent(self):
        """
        Get the pdb_details object that the fab belongs to.
        """
        return self.pdb

    def get_VH(self):
        """
        Get the heavy chain identifier for the fab.
        If this does not apply then 'NA' is returned
        """
        return self.VH

    def get_VL(self):
        """
        Get the light chain identifier for the fab.
        If this does not apply then 'NA' is returned
        """
        return self.VL
    
    def get_model(self):
        """
        Get the model number for the fab. 
        This is 0 for crystal structures.
        """
        return self.model
    
    def get_species(self):
        """
        Get the species of the Fab.
        """
        return self.species

    def get_common_species(self):
        """
        Get the common name species of the Fab.
        """
        return self.common_species
        
    def get_heavy_subclass(self):
        """
        Return the heavy chain subclass as defined by IMGT.
        As IMGT does not contain most of the data in our database this is 'unknown' for many fabs.
        """
        return self.heavy_subclass

    def get_light_subclass(self):
        """
        Return the light chain subclass as defined by IMGT.
        As IMGT does not contain most of the data in our database this is 'unknown' for many fabs.
        """
        return self.light_subclass
    
    def get_light_chain_type(self):
        """
        Return the light chain type (KAppa or Lambda) as defined by IMGT.
        As IMGT does not contain most of the data in our database this is 'unknown' for many fabs.
        """
        return self.light_chain_type

    def get_antigen(self):
        """
        Return the antigen_details object if an antigen is bound. Otherwise None.
        """
        return self.antigen

    def get_numbering(self,scheme=None):  # *  
        """
        Get the numbering for the heavy and light regions of the chain.
        """
        if scheme is None:
            scheme = self.pdb._database.get_numbering_scheme()
        if scheme != self.pdb._loaded_numbering_scheme:
            self._set_numbering(scheme=scheme)
        elif not hasattr(self, "numbering"):
            self._set_numbering(scheme=scheme)
        elif not self.numbering:
            self._set_numbering(scheme=scheme)        
            
        # this may return incorrectly if there are more than one regions on the chain...
        return self.numbering

    def get_sequence(self,scheme=None):  # *
        """
        Method to get the sequence of an fv. Note I put a chain break "/" to separate the two chains (or regions for scfv) 
        It will be VH and then VL
        If there is no VH or no VL the sequence will be represented by a gap character "-"
        """
        if scheme is None:
            scheme = self.pdb._database.get_numbering_scheme()
        if scheme != self.pdb._loaded_numbering_scheme:
            self._set_sequence(scheme=scheme)
        elif not hasattr(self, "sequence"):
            self._set_sequence(scheme=scheme)
        elif not self.sequence:
            self._set_sequence(scheme=scheme)

        return self.sequence
    
    def get_CDR_lengths(self, cdr=None, definition="chothia", scheme=None):
        """
        Count the number of residues in the CDRs. 
        
        The residues in the SEQRES SEQUENCE are given (not the structure which may have missing residues)
        
        @param cdr: A specific cdr to return the length of. e.g. Use CDRH3 or H3.
        @param defintion: The CDR definition to use. Choose one of "chothia","kabat" or "contact" or "imgt" or "north". "chothia" is used by default.
        @return: An integer if cdr is supplied, None if cdr is not present or recognised. Dictionary of CDR lengths otherwise.
        
        """
        definition = definition.lower()
        assert definition in ["chothia", "kabat", "contact", "imgt", "north", "wolfguy"], "Unrecognised CDR definition: %s" % definition

        if scheme is None:
            scheme = self.pdb._database.get_numbering_scheme()
        if scheme != self.pdb._loaded_numbering_scheme:
            self._set_CDR_lengths(scheme=scheme)
        elif not hasattr(self, "CDR_lengths"):
            self._set_CDR_lengths(scheme=scheme)
        elif not self.CDR_lengths:
            self._set_CDR_lengths(scheme=scheme)

        if type(cdr) is str and cdr:
            cdr=cdr.upper()
            if not cdr.startswith("CDR"):
                cdr = "CDR" + cdr.strip()
            try:
                return self.CDR_lengths[definition][cdr]
            except:
                AssertionError("Unrecognised CDR: %s"%cdr)
        else:
            return self.CDR_lengths[definition]

    def get_CDR_sequences(self, cdr=None, definition="chothia", scheme=None):
        """
        Get the annotated sequence of the CDRs

        The residues in the SEQRES SEQUENCE are given (not the structure which may have missing residues)

        @param cdr: A specific cdr to return the length of. e.g. Use CDRH3 or H3.
        @param defintion: The CDR definition to use. Choose one of "chothia","kabat" or "contact" or "imgt" or "north". "chothia" is used by default.        
        @return: An ordered list of (position, amino acid) tuples if cdr is supplied, None if cdr is not present or recognised. Dictionary of  ordered list of (position, amino acid) tuples otherwise
        """
        definition = definition.lower()
        assert definition in ["chothia", "kabat", "contact","imgt", "north","*"], "Unrecognised CDR definition: %s" % definition

        if scheme is None:
            scheme = self.pdb._database.get_numbering_scheme()
        if scheme != self.pdb._loaded_numbering_scheme:
            self._set_CDR_sequences(scheme=scheme)
        elif not hasattr(self, "CDR_sequences"):
            self._set_CDR_sequences(scheme=scheme)
        elif not self.CDR_sequences:
            self._set_CDR_sequences(scheme=scheme)

        if definition=="*": return self.CDR_sequences

        if type(cdr) is str and cdr:
            cdr=cdr.upper()
            if not cdr.startswith("CDR"):
                cdr = "CDR" + cdr.strip()
            try:
                return self.CDR_sequences[definition][cdr]
            except:
                raise AssertionError("CDR not recognised: %s"%cdr)
        else:
            return self.CDR_sequences[definition]

    def get_CDR_loops(self, cdr = None, definition = "chothia", scheme = None):
        """
        Get the joined (non-annotated) sequence of the CDRs
        """
        cdrs = self.get_CDR_sequences(cdr,definition,scheme)
        if cdr:
            return "".join([ y[1] for y in cdrs ])
        else:
            return dict([ (k, "".join([z[1] for z in v])) for k,v in list(cdrs.items()) ])

    def get_frameworks(self, definition = "chothia", scheme = None):
        """
        Get the joined (non-annotated) sequence of the framework regions
        """
        if hasattr(self, "fw"):
            return self.fw

        cdrs = self.get_CDR_sequences(definition = definition, scheme = scheme)
        numb = self.get_numbering(scheme=scheme)

        self.fw = {}

        if "H" in numb:
            heavynumb = set(numb["H"]) - (set([ y[0] for y in cdrs["CDRH1"]]) | set([z[0] for z in cdrs["CDRH2"]]) | set([x[0] for x in cdrs["CDRH3"]]) )
            self.fw["H"] = "".join([ numb["H"][n] for n in sorted(heavynumb, key = lambda h: (h[0], h[1])) ])
        if "L" in numb:
            lightnumb = set(numb["L"]) - (set([ y[0] for y in cdrs["CDRL1"]]) | set([z[0] for z in cdrs["CDRL2"]]) | set([x[0] for x in cdrs["CDRL3"]]) )
            self.fw["L"] = "".join([ numb["L"][n] for n in sorted(lightnumb, key = lambda h: (h[0], h[1])) ])

        return self.fw

    def get_CDR_clusters(self, cdr=None, definition="chothia", cutoff="0.5"):
        """
        Get the cdr cluster for a cdr or all of them.
        @param cdr: The cdr you wish to return. Leave as None to return all Choose from H1,H2,H3,L1,L2 or L3
        @param definition: The CDR characterisation you wish to use. Choose from chothia, kabat or contact
        @param cutoff: The UPGMA cutoff that should be used. This ensures no two structures in a cluster have an RMSD greater than this value. Choose from 0.5, 0.75,1.0 and 1.5

        The cluster returned will be the cluster id according to the sabdab clustering for the length.
        This will change each time more structures are added.
        Use get_same_cluster_structures to get structures with similar loops.
        """
        defs = ["kabat","chothia","contact"]
        cdrs = ["H1","H2","H3","L1","L2","L3"]
        cutoffs=["0.5","0.75","1.0","1.5"]
        assert definition in defs
        assert not cdr or cdr in cdrs
        assert cutoff in cutoffs
        if not hasattr(self, "CDR_clusters"):
            self.CDR_clusters={}
            self.other_structures_in_CDR_cluster={}
        if not self.pdb.clusteringpath:
            print("No clustering data found for the installed database", file=sys.stderr)
            return self.CDR_clusters
        try:
            if cdr:
                return self.CDR_clusters[cutoff][definition][cdr] # if loaded already
            else:
                return self.CDR_clusters[cutoff][definition]# if loaded already
        except KeyError:
            # descend the hierarchy
            if cutoff not in self.CDR_clusters:
               self.CDR_clusters[cutoff]={}
               self.other_structures_in_CDR_cluster[cutoff]={}
            if definition not in self.CDR_clusters[cutoff]:
                self.CDR_clusters[cutoff][definition]={}
                self.other_structures_in_CDR_cluster[cutoff][definition]={}
            for CDR in cdrs:
                if cdr and CDR!=cdr:continue # if they only want the one only populate the one.
                length=self.get_CDR_lengths(CDR, definition)
                cluster,others=None,[]
                if length:
                    filename="%s_%s%s_%d_CDR%s.pdb"%( self.id +( CDR, ))
                    clusteringfile = os.path.join(self.pdb.clusteringpath,"clusters_"+cutoff,definition,"CDR"+CDR,str(length)+"_clust.txt")
                    if os.path.exists(clusteringfile):
                        with open(clusteringfile) as openfile:
                            ordinal = 0
                            for line in openfile:
                                line = line.strip()
                                if "Cluster" in line:
                                    ordinal+=1
                                else:
                                    items = line.split(",")
                                    for item in items:
                                        if item==filename:
                                            cluster=ordinal
                                            break
                                    if cluster:
                                        for item in items:
                                            pdb, HL,model=item.split("_")[:3]
                                            if len(HL)==3:
                                                if HL.index("NA"): H,L=HL[0],"NA"
                                                else: H,L="NA",HL[-1]
                                            else:
                                                H,L=HL
                                            others.append((pdb,H,L,int(model)))
                                if cluster:break
                self.CDR_clusters[cutoff][definition][CDR]=cluster
                self.other_structures_in_CDR_cluster[cutoff][definition][CDR]=others
            if cdr:
                return self.CDR_clusters[cutoff][definition][cdr]
            else:
                return self.CDR_clusters[cutoff][definition]

    def get_same_cluster_structures(self,cdr=None, definition="chothia", cutoff="0.5"):
        """
        Get the identifiers of structures with the same CDR cluster as this fab.

        @param cdr: The cdr you wish to return. Leave as None to return all Choose from H1,H2,H3,L1,L2 or L3
        @param definition: The CDR characterisation you wish to use. Choose from chothia, kabat or contact
        @param cutoff: The UPGMA cutoff that should be used. This ensures no two structures in a cluster have an RMSD greater than this value. Choose from 0.5, 0.75,1.0 and 1.5

        The returned list (or dictionary of lists) will contain the structures in the database that are in the same cluster. 
        """
        if self.get_CDR_clusters(cdr, definition, cutoff):
            if cdr:
                return self.other_structures_in_CDR_cluster[cutoff][definition][cdr]
            else:
                return self.other_structures_in_CDR_cluster[cutoff][definition]
        elif cdr:
            return []
        else:
            return {}
        
    def get_affinity(self):
        """
        Get the affinity (Kd) of the antibody to its antigen. 
        All fabs in the pdb will have the same value if available (None otherwise)
                
        @return: The affinity (Kd) of the antibody to its antigen. Units are M. 
        @rtype: C{float} or None
        """
        return self.pdb.get_affinity()
    
    def get_orientation_angles(self, angle=""):
        """
        Get the orientation angles of the fv region.
        These are the angles between the VH-VL domain as described in:
        ABangle: Characterising the VH-VL orientation in antibodies. Dunbar J, Fuchs A, Shi J and Deane CM. PEDS. 2013
        
        @param angle: Optional argument for the angle returned. Choose one of "HL","HC1","HC2","LC1","LC2" or "dc".
        @return: Either a dictionary of the orientation angles or the specific value of the angle request.
        @rtype: C{float} or C{dict} of floats. Empty C{dict} if orientation angles cannot be found/calculated.
        """
        if not hasattr(self, "orientation_angles"):
            self._set_orientation_angles()
        if self.orientation_angles:
            if angle:
                return self.orientation_angles[angle]
            else:
                return self.orientation_angles
        else:
            return self.orientation_angles
        
    def get_imgt_details(self, chain_type=""):
        """
        Get the imgt details for the Fab.
        By default both the heavy and light imgt details dictionaries are returned.
        Use chain_type= H or L to retrieve for the chain.
        
        @param chain_type: Either "H" or "L" or "" for both.
        @return: A dictionary or tuple of Heavy Light dictionaries containing the imgt details for the chain.
        """
        if not hasattr(self, "heavy_imgt_details") or not hasattr(self, "light_imgt_details"):
            self._set_imgt_details()
        
        if chain_type.upper() == "H":
            return self.heavy_imgt_details
        elif chain_type.upper() == "L":
            return self.light_imgt_details
        elif not chain_type:
            return self.heavy_imgt_details, self.light_imgt_details
        else:
            raise Exception("Unrecognised chain type %s" % chain_type)

    def get_missing_residues(self, scheme=None):
        """
        Get those positions in the *variable region* that are missing in the structure.
        """

        if scheme is None:
            scheme = self.pdb._database.get_numbering_scheme()

        if scheme != self.pdb._loaded_numbering_scheme:
            self._set_missing_residues(scheme=scheme)
        elif not hasattr(self, "missing_residues"):
            self._set_missing_residues(scheme=scheme)

        return self.missing_residues

    def get_structure(self,scheme=None,definition=None):
        """
        Get the structure object for the fab. 
        """
        try:
            if self.VH == "NA":
                return self.pdb.get_structure(scheme=scheme,definition=definition)[int(self.model)][self.VL] # single chain
            elif self.VL == "NA": 
                return self.pdb.get_structure(scheme=scheme,definition=definition)[int(self.model)][self.VH] # single chain
            elif self.VH == self.VL:
                return self.pdb.get_structure(scheme=scheme,definition=definition)[int(self.model)][self.VH + self.VL.lower()] # scFv
            else: 
                return self.pdb.get_structure(scheme=scheme,definition=definition)[int(self.model)][self.VH+self.VL] # fab    
        except Exception as e: # lazy exception handling
            print("Problem loading structure for fab: %s %s %s %d"%self.id, file=sys.stderr)
            raise e
        
    # JD 13-03-15
    def get_identicals(self, threshold=0.99):
        """
        Get database entries that are identical (or within a sequence identity threshold to the fab)
        """
        threshold = float( threshold )
        assert threshold <= 1.0, "Sequence identity threshold must be <= 1.0"
        if not hasattr(self, "identicals"):
            self._set_identicals(threshold=threshold)
        elif "%.2f"%threshold != self._loaded_threshold:
            self._set_identicals(threshold=threshold)
        return self.identicals


    # JL 12-03-15
    # JD 13-03-15 mod
    def get_unbound_forms(self,threshold = 0.99):
        """
        Retrieve the other unbound forms of this fv from the database
        Currently only available for paired VH and VL

        @param threshold: The sequence identity theshold that should be used
        @return: A dictionary of fab objects that are the same and unbound.
        """            
        assert self.VH != "NA" and self.VL != "NA", "Currently only available for paired VH and VL"
        threshold = float( threshold )
        assert threshold <= 1.0, "Sequence identity threshold must be <= 1.0"
        if not hasattr(self, "unbound_forms"):
            self._set_unbound_forms(threshold=threshold)
        elif "%.2f"%threshold != self._loaded_threshold:
            self._set_unbound_forms(threshold=threshold)
        return self.unbound_forms

    # JL 12-03-15
    # JD 13-03-15 mod
    def get_bound_forms(self,threshold = 0.99):
        """
        Retrieve the other bound forms of this fv from the database.
        Currently only available for paired VH and VL
        
        @param threshold: The sequence identity theshold that should be used
        @return: A dictionary of fab objects that are the same and bound.
        """            
        assert self.VH != "NA" and self.VL != "NA", "Currently only available for paired VH and VL"
        threshold = float( threshold )
        assert threshold <= 1.0, "Sequence identity threshold must be <= 1.0"
        if not hasattr(self, "bound_forms"):
            self._set_bound_forms(threshold=threshold)
        elif "%.2f"%threshold != self._loaded_threshold:
            self._set_bound_forms(threshold=threshold)
        return self.bound_forms

    ##JL
    ##def primeNumberGenerator(self):
        ## return 571


    def predict_north_canonical( self, cdr, scheme=None ):
        """
        Predict the canonical ID, representative's ID and the sequence.
        """
        sequence = self.get_CDR_sequences( cdr, definition="north", scheme=scheme )
        return assign_canonical( "".join( [ _[1] for _ in sequence]), cdr )


    def is_scfv(self):
        """
        Is the fab a single chain fv?
        """
        return self.scfv

    def is_complex(self):
        """
        Is the fab in complex with an antigen?
        """
        return self.complex
    
    def is_completefab(self):
        """
        Is the fab a completefab (i.e. is it paired?)
        """
        return self.completefab

    def has_constant(self):
        """
        Does the fab have constant domains? 
        """
        if "C" in self.get_imgt_details("H") or "C" in self.get_imgt_details("L"):
            return True
        else:
            return False

       
    def get_summary_line(self):
        """
        Get a tab-separated line for fields for the fab entry.
        """
        try:
            return self.pdb.get_summary_lines()[(self.pdb.id, self.VH, self.VL, str(self.model))]
        except KeyError:
            return ""
        
    
class Antigen_details(object):
    '''
    Class to hold the details of an antigen
    '''
    def __init__(self, chain=None):
        """
        Initialise the antigen with the chain that it is on (or is)
        
        """
        self._set_chains(chain)
        self.sequence = ""
        self.sequences={}
        
    def __str__(self):
        if self.chains:
            return "antigen chain = %s" % " | ".join(self.chains)
        else:
            return ""

    def __repr__(self):
        if self.chains:
            return "\nantigen chain = %s\n" % " | ".join(self.chains)
        else:
            return "None"
    
    def _set_chains(self, chain=None):
        """
        Set the chain that the antigen is on.
        If the antigen is multi-chained then all details will be separated by " | "
        """
        self.chains=[c.strip() for c in chain.split("|")]
#        self.chain = chain
             
    def _set_antigen_type(self, agtype=None):
        """
        Protein,peptide, hapten, nucleic acid, carbohydrate
        """
        if self.is_multichain():
            if "protein" in agtype:
                self.agtype="protein"
            elif "peptide" in agtype:
                self.agtype="peptide"
            elif "nucleic-acid" in agtype:
                self.agtype="nucleic-acid"
            elif "hapten" in agtype:
                self.agtype="hapten"
            else:
                self.agtype="unknown"
        else:
            self.agtype = agtype
        
    def _set_antigen_het_name(self, hetname=None):
        """
        Set the het name of the antigen.
        This will be NA if the antigen is formed of ATOM records.
        Applies mainly for hapten antigens.
        """
        self.het_name = hetname

    def _set_antigen_name(self, name=None):
        """
        Set the name of the antigen.
        This has either been parsed from the chemical component database or from the pdb header depending on the antigen type.
        @change: uniquify the names but preserve ordering
        """
        names = []
        for c in name.split("|"):
            if c.strip() not in names:
                names.append(c.strip())
        self.name = "; ".join( names )

    def _set_pdb(self, pdb):
        """
        Set the pdb that the antigen comes from
        """
        self.pdb = pdb  # *

    def _set_fab(self, fab):
        """
        Set the fab that the antigen is bound to
        """
        self.fab = fab  # *
        
    def _set_species(self, species=None):
        if species:
            s=list(map(str.strip,species.split("|")))
            if len(set(s)) > 1:
                self.species =  " ; ".join(s)
            else:
                self.species = list(set(s))[0].upper()  ####
        else:  # this is found from pdb taxonomy. -if it is empty then it is a synthetic construct.
            self.species = "SYNTHETIC CONSTRUCT"

    def get_chain(self):
        """
        Tells which chain the antigen is on.
        """
        if self.is_multichain():
            raise Exception("Use get_chains to find the chain ids of multi chain antigens") 
        return self.chains[0]

    def get_chains(self):
        """
        Tells which chain the antigen is on.
        """
        return self.chains

    def get_antigen_type(self):
        """
        Return the type of antigen.
        """
        return self.agtype

    def get_antigen_name(self):
        """
        Get the chemical name of the antigen.
        """
        return self.name

    def get_species(self):
        return self.species

    def get_antigen_het_name(self):
        """
        Get the hetname of the antigen if it was from a hetatm
        """
        return self.het_name

    def get_fab(self):
        """
        Get the fab that the antigen is bound to.
        """
        return self.fab
    
    def get_sequences(self):
        """
        Get the sequences of the antigen chains if it is a peptide or protein.
        """
        if self.get_antigen_type() in ["peptide", "protein"]:
            if not self.sequences:
                for chain in self.get_chains():
                    self.sequences[chain] = self.pdb.get_raw_sequences()[chain]
        else:
            print("Cannot return non peptide/protein antigen sequences yet.")
        return self.sequences

    
    def get_sequence(self):
        """
        Get the sequence of the antigen if it is a peptide or protein.
        If the antigen is multi-chained, then this function will return the sequence string of all the chains separated by "/"
        These chains appear in order of their length and then their hash.
        This should ensure if identical antigens appear in two complexes, then the sequences will appear the same. 
        """
        if self.get_antigen_type() in ["peptide", "protein"]:
            sequences = self.get_sequences()
            return "/".join( [ sequences[chain] for chain in sorted( sequences, key=lambda x: ( len( sequences[x]), sequences[x] ) )  ] )
        else:
            print("Cannot return non peptide/protein antigen sequences yet.")
        return self.sequence

        
    def is_multichain(self):
        if len(self.get_chains())>1:
            return True
        else:
            return False
