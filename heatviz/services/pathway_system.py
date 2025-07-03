# import os
# import json
# import re
# from collections import defaultdict
# import glob

# class PathwayDatabase:
#     def __init__(self, pathway_data_dir="pathway_data"):
#         """
#         Universal pathway database that handles multiple pathway library formats
        
#         Args:
#             pathway_data_dir (str): Directory containing pathway txt files
#         """
#         self.pathway_data_dir = pathway_data_dir
#         self.pathways = {}  # pathway_name -> pathway_info
#         self.genes_to_pathways = defaultdict(list)  # gene -> list of pathways
        
#         # Enhanced functional categories for all pathway types
#         self.categories = {
#             "immune": [
#                 "immune", "immunoglobulin", "interferon", "cytokine", "inflammation", 
#                 "tcell", "bcell", "macrophage", "dendritic", "lymphocyte", "antiviral", 
#                 "innate", "adaptive", "autoimmune", "allergy", "th1", "th2", "treg", 
#                 "nk_cell", "neutrophil", "eosinophil", "antibody", "complement", 
#                 "toll", "interleukin", "tumor_necrosis", "mhc", "antigen", "response"
#             ],
#             "cancer": [
#                 "cancer", "tumor", "oncogene", "suppressor", "metastasis", "apoptosis", 
#                 "cell_cycle", "proliferation", "angiogenesis", "invasion", "malignant",
#                 "carcinoma", "sarcoma", "leukemia", "lymphoma", "neoplasm", "growth",
#                 "dna_repair", "p53", "rb", "myc", "brca", "oncogenesis", "mitosis"
#             ],
#             "metabolism": [
#                 "metabolism", "metabolic", "glucose", "lipid", "fatty_acid", "glycolysis", 
#                 "gluconeogenesis", "oxidative", "mitochondria", "insulin", "diabetes", 
#                 "cholesterol", "amino_acid", "energy", "atp", "respiratory", "citric",
#                 "pentose", "lipogenesis", "lipolysis", "ketogenesis", "biosynthesis"
#             ],
#             "development": [
#                 "development", "developmental", "embryonic", "differentiation", "neural", 
#                 "mesoderm", "endoderm", "ectoderm", "organogenesis", "morphogenesis", 
#                 "patterning", "stem_cell", "pluripotency", "regeneration", "embryo",
#                 "fetal", "tissue", "organ", "limb", "brain", "heart", "kidney"
#             ],
#             "signaling": [
#                 "signaling", "signal", "transduction", "mapk", "pi3k", "akt", "mtor",
#                 "wnt", "notch", "hedgehog", "tgf", "bmp", "jak", "stat", "calcium",
#                 "camp", "cgmp", "protein_kinase", "phosphatase", "receptor", "cascade"
#             ],
#             "stress": [
#                 "stress", "oxidative", "endoplasmic", "reticulum", "hypoxia", "heat_shock", 
#                 "dna_damage", "unfolded", "protein", "detoxification", "antioxidant", 
#                 "autophagy", "senescence", "aging", "response", "repair", "damage"
#             ],
#             "transport": [
#                 "transport", "transmembrane", "channel", "transporter", "pump", "carrier",
#                 "ion", "sodium", "potassium", "calcium", "chloride", "membrane", 
#                 "vesicle", "endocytosis", "exocytosis", "trafficking", "localization"
#             ],
#             "transcription": [
#                 "transcription", "transcriptional", "rna", "polymerase", "promoter", 
#                 "enhancer", "silencer", "chromatin", "histone", "methylation", 
#                 "acetylation", "epigenetic", "gene_expression", "regulation"
#             ],
#             "enzyme": [
#                 "enzyme", "catalytic", "kinase", "phosphatase", "transferase", 
#                 "hydrolase", "ligase", "lyase", "isomerase", "oxidoreductase",
#                 "peptidase", "protease", "dehydrogenase", "synthase", "reductase"
#             ],
#             "structural": [
#                 "structural", "cytoskeleton", "actin", "tubulin", "intermediate", 
#                 "filament", "microtubule", "microfilament", "cell_junction", 
#                 "adhesion", "extracellular_matrix", "collagen", "elastin"
#             ]
#         }
        
#         self.load_database()

#     def detect_pathway_format(self, filename, sample_lines):
#         """
#         Detect the format of a pathway file based on filename and content
        
#         Args:
#             filename (str): Name of the file
#             sample_lines (list): First few lines of the file
            
#         Returns:
#             str: Format type
#         """
#         filename_lower = filename.lower()
        
#         if filename_lower.startswith("go_"):
#             return "gene_ontology"
#         elif filename_lower.startswith("chea_") or "chea" in filename_lower:
#             return "transcription_factor"
#         elif filename_lower.startswith("kegg_"):
#             return "kegg"
#         elif filename_lower.startswith("reactome_"):
#             return "reactome"
#         elif filename_lower.startswith("msigdb_"):
#             return "msigdb"
#         elif filename_lower.startswith("wikipathway_"):
#             return "wikipathway"
#         elif filename_lower.startswith("biocarta_"):
#             return "biocarta"
#         elif filename_lower.startswith("bioplanet_"):
#             return "bioplanet"
#         elif "gtex" in filename_lower:
#             return "gtex"
#         elif "atlas" in filename_lower:
#             return "atlas"
#         elif "disease" in filename_lower or "omim" in filename_lower:
#             return "disease"
#         elif "phenotype" in filename_lower:
#             return "phenotype"
#         elif "trrust" in filename_lower:
#             return "transcription_factor"
#         elif "encode" in filename_lower:
#             return "transcription_factor"
#         elif "dsigdb" in filename_lower:
#             return "drug"
#         else:
#             # Try to detect from content
#             if sample_lines:
#                 first_line = sample_lines[0]
#                 if "\t" in first_line:
#                     parts = first_line.split("\t")
#                     if len(parts) >= 5 and any(x.isdigit() for x in parts[1:3]):
#                         return "transcription_factor"  # Likely ChEA format
#                     elif "GO:" in first_line:
#                         return "gene_ontology"
#                     else:
#                         return "generic_gmt"  # GMT format
            
#             return "generic_gmt"

#     def parse_gene_ontology(self, file_path, library_name):
#         """Parse Gene Ontology format files"""
#         pathways_count = 0
        
#         with open(file_path, 'r', encoding='utf-8') as f:
#             for line_num, line in enumerate(f, 1):
#                 line = line.strip()
#                 if not line:
#                     continue
                
#                 try:
#                     # GO format: GO_Library_Function_Description_Gene1_Gene2_...
#                     parts = line.split('_')
#                     if len(parts) < 4:
#                         continue
                    
#                     # Extract GO function description (everything between library and genes)
#                     library_prefix = f"{parts[0]}_{parts[1]}_{parts[2]}"  # GO_Molecular_Function_2023
                    
#                     # Find where genes start (usually after GO term in parentheses)
#                     description_parts = []
#                     genes = []
#                     in_description = True
                    
#                     for part in parts[3:]:
#                         if in_description and ("GO:" in part or len(genes) > 0):
#                             in_description = False
                        
#                         if in_description:
#                             description_parts.append(part)
#                         else:
#                             if part and part not in description_parts:
#                                 genes.append(part)
                    
#                     if not genes:
#                         # If no genes found, treat last parts as genes
#                         genes = parts[-2:]  # Take last 2 as genes
#                         description_parts = parts[3:-2]
                    
#                     function_description = " ".join(description_parts)
#                     pathway_name = f"{library_prefix}_{function_description}"
                    
#                     if genes:
#                         self.pathways[pathway_name] = {
#                             'genes': genes,
#                             'description': function_description,
#                             'library': library_name,
#                             'type': 'gene_ontology',
#                             'gene_count': len(genes)
#                         }
                        
#                         # Index genes
#                         for gene in genes:
#                             self.genes_to_pathways[gene].append(pathway_name)
                        
#                         pathways_count += 1
                        
#                 except Exception as e:
#                     continue
        
#         return pathways_count

#     def parse_transcription_factor(self, file_path, library_name):
#         """Parse transcription factor format (ChEA, TRRUST, ENCODE)"""
#         pathways_count = 0
        
#         with open(file_path, 'r', encoding='utf-8') as f:
#             for line_num, line in enumerate(f, 1):
#                 line = line.strip()
#                 if not line:
#                     continue
                
#                 try:
#                     parts = line.split('\t')
#                     if len(parts) < 6:
#                         continue
                    
#                     tf_name = parts[0]
#                     pmid = parts[1] if len(parts) > 1 else ""
#                     method = parts[2] if len(parts) > 2 else ""
#                     cell_type = parts[3] if len(parts) > 3 else ""
#                     species = parts[4] if len(parts) > 4 else ""
#                     genes = [gene.strip() for gene in parts[5:] if gene.strip()]
                    
#                     if not genes:
#                         continue
                    
#                     pathway_name = f"{library_name}_{tf_name}_{cell_type}_{species}"
                    
#                     self.pathways[pathway_name] = {
#                         'genes': genes,
#                         'tf': tf_name,
#                         'cell_type': cell_type,
#                         'species': species,
#                         'method': method,
#                         'pmid': pmid,
#                         'library': library_name,
#                         'type': 'transcription_factor',
#                         'gene_count': len(genes)
#                     }
                    
#                     # Index genes
#                     for gene in genes:
#                         self.genes_to_pathways[gene].append(pathway_name)
                    
#                     pathways_count += 1
                    
#                 except Exception as e:
#                     continue
        
#         return pathways_count

#     def parse_generic_gmt(self, file_path, library_name):
#         """Parse generic GMT format (KEGG, Reactome, MSigDB, etc.)"""
#         pathways_count = 0
        
#         with open(file_path, 'r', encoding='utf-8') as f:
#             for line_num, line in enumerate(f, 1):
#                 line = line.strip()
#                 if not line:
#                     continue
                
#                 try:
#                     parts = line.split('\t')
#                     if len(parts) < 3:
#                         continue
                    
#                     pathway_name = parts[0]
#                     description = parts[1] if parts[1] else pathway_name
#                     genes = [gene.strip() for gene in parts[2:] if gene.strip()]
                    
#                     if not genes:
#                         continue
                    
#                     full_pathway_name = f"{library_name}_{pathway_name}"
                    
#                     self.pathways[full_pathway_name] = {
#                         'genes': genes,
#                         'description': description,
#                         'library': library_name,
#                         'type': 'pathway',
#                         'gene_count': len(genes)
#                     }
                    
#                     # Index genes
#                     for gene in genes:
#                         self.genes_to_pathways[gene].append(full_pathway_name)
                    
#                     pathways_count += 1
                    
#                 except Exception as e:
#                     continue
        
#         return pathways_count

#     def load_database(self):
#         """Load all pathway txt files using appropriate parsers"""
#         print(f"📚 Loading universal pathway database from {self.pathway_data_dir}/")
        
#         if not os.path.exists(self.pathway_data_dir):
#             print(f"❌ Directory {self.pathway_data_dir} not found.")
#             return
        
#         total_pathways = 0
#         txt_files = glob.glob(os.path.join(self.pathway_data_dir, "*.txt"))
        
#         for file_path in txt_files:
#             library_name = os.path.basename(file_path).replace('.txt', '')
            
#             # Read first few lines to detect format
#             try:
#                 with open(file_path, 'r', encoding='utf-8') as f:
#                     sample_lines = [f.readline().strip() for _ in range(3)]
#                     sample_lines = [line for line in sample_lines if line]
#             except:
#                 sample_lines = []
            
#             # Detect format
#             format_type = self.detect_pathway_format(library_name, sample_lines)
            
#             # Parse based on format
#             if format_type == "gene_ontology":
#                 pathways_loaded = self.parse_gene_ontology(file_path, library_name)
#             elif format_type == "transcription_factor":
#                 pathways_loaded = self.parse_transcription_factor(file_path, library_name)
#             else:
#                 pathways_loaded = self.parse_generic_gmt(file_path, library_name)
            
#             total_pathways += pathways_loaded
#             print(f"  ✅ Loaded {pathways_loaded} pathways from {library_name} ({format_type})")
        
#         print(f"📊 Total pathways loaded: {total_pathways}")
#         print(f"📊 Total unique genes indexed: {len(self.genes_to_pathways)}")

#     def search_pathways_by_category(self, category):
#         """Search pathways by functional category"""
#         print(f"🔍 Searching for '{category}' pathways...")
        
#         category_lower = category.lower()
#         matches = []
#         seen_pathways = set()
        
#         # Get keywords for this category
#         keywords = self.categories.get(category_lower, [category_lower])
#         print(f"Using keywords: {keywords[:5]}...")
        
#         # Search in pathway names and descriptions
#         for pathway_name, pathway_data in self.pathways.items():
#             if pathway_name in seen_pathways:
#                 continue
            
#             # Create searchable text
#             search_text_parts = [
#                 pathway_name,
#                 pathway_data.get('description', ''),
#                 pathway_data.get('tf', ''),
#                 pathway_data.get('cell_type', ''),
#                 pathway_data.get('library', '')
#             ]
#             search_text = ' '.join(str(part) for part in search_text_parts).lower()
            
#             # Check if any keyword matches
#             for keyword in keywords:
#                 if keyword in search_text:
#                     matches.append({
#                         'pathway_name': pathway_name,
#                         'description': pathway_data.get('description', pathway_name),
#                         'library': pathway_data['library'],
#                         'type': pathway_data['type'],
#                         'gene_count': pathway_data['gene_count'],
#                         'match_reason': f"Contains '{keyword}'",
#                         'confidence': 'high' if keyword == category_lower else 'medium'
#                     })
#                     seen_pathways.add(pathway_name)
#                     break
        
#         # Sort by confidence and gene count
#         matches.sort(key=lambda x: (
#             x['confidence'] != 'high',
#             -x['gene_count']
#         ))
        
#         print(f"Found {len(matches)} matching pathways")
#         return matches[:30]

#     def get_pathway_genes(self, pathway_query, action_type="pathway_filter"):
#         """Get genes for a specific pathway"""
#         print(f"🔍 Looking for pathway: '{pathway_query}'")
        
#         # Try exact match first
#         if pathway_query in self.pathways:
#             genes = self.pathways[pathway_query]['genes']
#             print(f"Found exact match: {len(genes)} genes")
#             return genes
        
#         # Try partial match in pathway names
#         pathway_lower = pathway_query.lower()
#         for pathway_name, pathway_data in self.pathways.items():
#             if pathway_lower in pathway_name.lower():
#                 genes = pathway_data['genes']
#                 print(f"Found partial match '{pathway_name}': {len(genes)} genes")
#                 return genes
        
#         # Try searching in descriptions
#         for pathway_name, pathway_data in self.pathways.items():
#             description = pathway_data.get('description', '').lower()
#             if pathway_lower in description:
#                 genes = pathway_data['genes']
#                 print(f"Found in description '{pathway_name}': {len(genes)} genes")
#                 return genes
        
#         print(f"No pathway found for: '{pathway_query}'")
#         return []

#     def get_functional_genes(self, function_name):
#         """Get genes associated with a biological function"""
#         matching_pathways = self.search_pathways_by_category(function_name)
        
#         all_genes = []
#         for pathway in matching_pathways:
#             pathway_genes = self.get_pathway_genes(pathway['pathway_name'])
#             all_genes.extend(pathway_genes)
        
#         unique_genes = list(set(all_genes))
#         print(f"Found {len(unique_genes)} unique genes for function '{function_name}'")
#         return unique_genes

# # Global database instance
# pathway_db = None

# def initialize_pathway_database(pathway_data_dir="pathway_data"):
#     """Initialize the global pathway database"""
#     global pathway_db
#     if pathway_db is None:
#         pathway_db = PathwayDatabase(pathway_data_dir)
#     return pathway_db

# def get_pathway_genes(pathway_name, action_type="pathway_filter"):
#     """Get genes for a pathway"""
#     global pathway_db
#     if pathway_db is None:
#         pathway_db = initialize_pathway_database()
#     return pathway_db.get_pathway_genes(pathway_name, action_type)

# def search_pathways_by_category(category):
#     """Search pathways by category"""
#     global pathway_db
#     if pathway_db is None:
#         pathway_db = initialize_pathway_database()
#     return pathway_db.search_pathways_by_category(category)

# def get_functional_genes(function_name):
#     """Get genes for a biological function"""
#     global pathway_db
#     if pathway_db is None:
#         pathway_db = initialize_pathway_database()
#     return pathway_db.get_functional_genes(function_name)

# def validate_pathway_command(action, target, value):
#     """Validate pathway commands"""
#     global pathway_db
#     if pathway_db is None:
#         pathway_db = initialize_pathway_database()
    
#     if action == "pathway_filter":
#         genes = pathway_db.get_pathway_genes(value)
#         if not genes:
#             return False, f"Pathway '{value}' not found. Try 'list {value} pathways' to see available options."
    
#     elif action == "pathway_search":
#         results = pathway_db.search_pathways_by_category(value)
#         if not results:
#             available_categories = list(pathway_db.categories.keys())
#             return False, f"No pathways found for '{value}'. Available categories: {', '.join(available_categories[:5])}, etc."
    
#     elif action == "functional_filter":
#         genes = pathway_db.get_functional_genes(value)
#         if not genes:
#             return False, f"No genes found for function '{value}'. Try searching for related pathways first."
    
#     return True, ""

# # Example usage
# if __name__ == "__main__":
#     db = PathwayDatabase("pathway_data")
    
#     print("\n🧪 Testing searches:")
    
#     # Test immune search
#     immune_results = db.search_pathways_by_category("immune")
#     print(f"Immune pathways: {len(immune_results)}")
#     for i, result in enumerate(immune_results[:3]):
#         print(f"  {i+1}. {result['description'][:60]}... ({result['gene_count']} genes)")
    
#     # Test histone search
#     histone_results = db.search_pathways_by_category("histone")
#     print(f"\nHistone pathways: {len(histone_results)}")
#     for i, result in enumerate(histone_results[:3]):
#         print(f"  {i+1}. {result['description'][:60]}... ({result['gene_count']} genes)")


import os
import json
import re
from collections import defaultdict
import glob

class PathwayDatabase:
    def __init__(self, pathway_data_dir="pathway_data"):
        """
        Focused pathway database for biological pathways only
        
        Args:
            pathway_data_dir (str): Directory containing pathway txt files
        """
        self.pathway_data_dir = pathway_data_dir
        self.pathways = {}  # pathway_name -> pathway_info
        self.genes_to_pathways = defaultdict(list)  # gene -> list of pathways
        
        # Focused pathway categories (removed transcription & ontology categories)
        self.categories = {
            "cancer": [
                "cancer", "tumor", "oncogene", "suppressor", "metastasis", "apoptosis", 
                "cell_cycle", "proliferation", "angiogenesis", "invasion", "malignant",
                "carcinoma", "sarcoma", "leukemia", "lymphoma", "neoplasm", "growth",
                "dna_repair", "p53", "rb", "myc", "brca", "oncogenesis", "mitosis"
            ],
            "immune": [
                "immune", "immunoglobulin", "interferon", "cytokine", "inflammation", 
                "tcell", "bcell", "macrophage", "dendritic", "lymphocyte", "antiviral", 
                "innate", "adaptive", "autoimmune", "allergy", "th1", "th2", "treg", 
                "nk_cell", "neutrophil", "eosinophil", "antibody", "complement", 
                "toll", "interleukin", "tumor_necrosis", "mhc", "antigen", "response"
            ],
            "metabolism": [
                "metabolism", "metabolic", "glucose", "lipid", "fatty_acid", "glycolysis", 
                "gluconeogenesis", "oxidative", "mitochondria", "insulin", "diabetes", 
                "cholesterol", "amino_acid", "energy", "atp", "respiratory", "citric",
                "pentose", "lipogenesis", "lipolysis", "ketogenesis", "biosynthesis"
            ],
            "signaling": [
                "signaling", "signal", "transduction", "mapk", "pi3k", "akt", "mtor",
                "wnt", "notch", "hedgehog", "tgf", "bmp", "jak", "stat", "calcium",
                "camp", "cgmp", "protein_kinase", "phosphatase", "receptor", "cascade"
            ],
            "development": [
                "development", "developmental", "embryonic", "differentiation", "neural", 
                "mesoderm", "endoderm", "ectoderm", "organogenesis", "morphogenesis", 
                "patterning", "stem_cell", "pluripotency", "regeneration", "embryo",
                "fetal", "tissue", "organ", "limb", "brain", "heart", "kidney"
            ],
            "stress": [
                "stress", "oxidative", "endoplasmic", "reticulum", "hypoxia", "heat_shock", 
                "dna_damage", "unfolded", "protein", "detoxification", "antioxidant", 
                "autophagy", "senescence", "aging", "response", "repair", "damage"
            ],
            "transport": [
                "transport", "transmembrane", "channel", "transporter", "pump", "carrier",
                "ion", "sodium", "potassium", "calcium", "chloride", "membrane", 
                "vesicle", "endocytosis", "exocytosis", "trafficking", "localization"
            ],
            "disease": [
                "disease", "disorder", "syndrome", "pathology", "dysfunction", "infection",
                "viral", "bacterial", "fungal", "cardiovascular", "neurological", 
                "psychiatric", "autoimmune", "inflammatory", "genetic", "hereditary"
            ]
        }
        
        # Only pathway libraries to include
        self.supported_libraries = {
            "kegg_2021_human",
            "wikipathways_2023_human", 
            "reactome_2022",
            "biocarta_2016",
            "msigdb_hallmark_2020",
            "nci_nature_2016",
            "panther_2016"
        }
        
        self.load_database()

    def is_pathway_library(self, filename):
        """Check if file is from a supported pathway library"""
        filename_lower = filename.lower()
        
        # Your exact pathway libraries (lowercase to match filename.lower())
        supported_files = [
            "biocarta_2016",
            "kegg_2021_human", 
            "msigdb_hallmark_2020",
            "reactome_2022",
            "wikipathways_2023_human"  # Note: "wikipathways" not "wikipathway"
        ]
        
        # Check exact matches first
        for supported in supported_files:
            if supported in filename_lower:
                return True
        
        # Fallback: general pathway indicators
        pathway_indicators = [
            "kegg", "wikipathway", "reactome", "biocarta", "msigdb_hallmark"
        ]
        
        # Skip transcription and ontology libraries  
        skip_indicators = [
            "go_", "chea_", "encode", "trrust", "gtex", "atlas",
            "histone", "chip", "transcription", "tf_", "mirna"
        ]
        
        # Check if should be skipped
        for skip in skip_indicators:
            if skip in filename_lower:
                return False
        
        # Check if is pathway library
        for pathway in pathway_indicators:
            if pathway in filename_lower:
                return True
        
        return False

    def detect_pathway_format(self, filename, sample_lines):
        """
        Detect pathway format - simplified to only handle pathway formats
        """
        filename_lower = filename.lower()
        
        if filename_lower.startswith("kegg_"):
            return "kegg"
        elif filename_lower.startswith("reactome_"):
            return "reactome"
        elif filename_lower.startswith("msigdb_"):
            return "msigdb"
        elif filename_lower.startswith("wikipathway_"):
            return "wikipathway"
        elif filename_lower.startswith("biocarta_"):
            return "biocarta"
        elif filename_lower.startswith("nci_"):
            return "nci"
        elif filename_lower.startswith("panther_"):
            return "panther"
        else:
            return "generic_pathway"

    def parse_pathway_gmt(self, file_path, library_name):
        """Parse pathway GMT format files"""
        pathways_count = 0
        
        with open(file_path, 'r', encoding='utf-8') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line:
                    continue
                
                try:
                    parts = line.split('\t')
                    if len(parts) < 3:
                        continue
                    
                    pathway_name = parts[0]
                    description = parts[1] if parts[1] else pathway_name
                    genes = [gene.strip() for gene in parts[2:] if gene.strip()]
                    
                    if not genes:
                        continue
                    
                    # Clean pathway name
                    clean_pathway_name = pathway_name.replace(" ", "_").replace("-", "_")
                    full_pathway_name = f"{library_name}_{clean_pathway_name}"
                    
                    self.pathways[full_pathway_name] = {
                        'genes': genes,
                        'description': description,
                        'original_name': pathway_name,
                        'library': library_name,
                        'type': 'pathway',
                        'gene_count': len(genes)
                    }
                    
                    # Index genes
                    for gene in genes:
                        self.genes_to_pathways[gene].append(full_pathway_name)
                    
                    pathways_count += 1
                    
                except Exception as e:
                    print(f"Warning: Error parsing line {line_num} in {library_name}: {e}")
                    continue
        
        return pathways_count

    def load_database(self):
        """Load only pathway txt files"""
        print(f"📚 Loading pathway-focused database from {self.pathway_data_dir}/")
        
        if not os.path.exists(self.pathway_data_dir):
            print(f"❌ Directory {self.pathway_data_dir} not found.")
            return
        
        total_pathways = 0
        total_skipped = 0
        txt_files = glob.glob(os.path.join(self.pathway_data_dir, "*.txt"))
        
        for file_path in txt_files:
            library_name = os.path.basename(file_path).replace('.txt', '')
            
            # Only process pathway libraries
            if not self.is_pathway_library(library_name):
                print(f"  ⏭️  Skipping {library_name} (not a pathway library)")
                total_skipped += 1
                continue
            
            # Read first few lines to detect format
            try:
                with open(file_path, 'r', encoding='utf-8') as f:
                    sample_lines = [f.readline().strip() for _ in range(3)]
                    sample_lines = [line for line in sample_lines if line]
            except:
                sample_lines = []
            
            # Detect format
            format_type = self.detect_pathway_format(library_name, sample_lines)
            
            # Parse pathway file
            pathways_loaded = self.parse_pathway_gmt(file_path, library_name)
            
            total_pathways += pathways_loaded
            print(f"  ✅ Loaded {pathways_loaded} pathways from {library_name} ({format_type})")
        
        print(f"📊 Total pathways loaded: {total_pathways}")
        print(f"📊 Total files skipped: {total_skipped}")
        print(f"📊 Total unique genes indexed: {len(self.genes_to_pathways)}")
        print(f"📊 Pathway libraries loaded: {len([f for f in txt_files if self.is_pathway_library(os.path.basename(f).replace('.txt', ''))])}")

    def search_pathways_by_category(self, category):
        """Search pathways by functional category"""
        print(f"🔍 Searching for '{category}' pathways...")
        
        category_lower = category.lower()
        matches = []
        seen_pathways = set()
        
        # Get keywords for this category
        keywords = self.categories.get(category_lower, [category_lower])
        
        # Add common variations of the search term
        if category_lower not in keywords:
            keywords.append(category_lower)
        
        print(f"Using keywords: {keywords[:5]}{'...' if len(keywords) > 5 else ''}")
        
        # Search in pathway names and descriptions
        for pathway_name, pathway_data in self.pathways.items():
            if pathway_name in seen_pathways:
                continue
            
            # Create searchable text from pathway info
            search_text_parts = [
                pathway_data.get('original_name', ''),
                pathway_data.get('description', ''),
                pathway_name
            ]
            search_text = ' '.join(str(part) for part in search_text_parts).lower()
            
            # Check if any keyword matches
            match_found = False
            matched_keyword = ""
            for keyword in keywords:
                if keyword in search_text:
                    match_found = True
                    matched_keyword = keyword
                    break
            
            if match_found:
                # Calculate confidence based on keyword match
                confidence = 'high' if matched_keyword == category_lower else 'medium'
                
                matches.append({
                    'pathway_name': pathway_name,
                    'original_name': pathway_data.get('original_name', pathway_name),
                    'description': pathway_data.get('description', pathway_name),
                    'library': pathway_data['library'],
                    'type': pathway_data['type'],
                    'gene_count': pathway_data['gene_count'],
                    'match_reason': f"Contains '{matched_keyword}'",
                    'confidence': confidence
                })
                seen_pathways.add(pathway_name)
        
        # Sort by confidence, then by gene count (larger pathways first)
        matches.sort(key=lambda x: (
            x['confidence'] != 'high',  # High confidence first
            x['confidence'] != 'medium',  # Then medium confidence
            -x['gene_count']  # Then by gene count (descending)
        ))
        
        print(f"Found {len(matches)} matching pathways")
        return matches[:20]  # Return top 20 matches

    def get_pathway_genes(self, pathway_query, action_type="pathway_filter"):
        """Get genes for a specific pathway"""
        print(f"🔍 Looking for pathway: '{pathway_query}'")
        
        # Try exact match first
        if pathway_query in self.pathways:
            genes = self.pathways[pathway_query]['genes']
            print(f"Found exact match: {len(genes)} genes")
            return genes
        
        # Try partial match in pathway names
        pathway_lower = pathway_query.lower()
        best_match = None
        best_score = 0
        
        for pathway_name, pathway_data in self.pathways.items():
            # Check pathway name
            if pathway_lower in pathway_name.lower():
                score = len(pathway_lower) / len(pathway_name)
                if score > best_score:
                    best_score = score
                    best_match = (pathway_name, pathway_data)
            
            # Check original name
            original_name = pathway_data.get('original_name', '').lower()
            if original_name and pathway_lower in original_name:
                score = len(pathway_lower) / len(original_name)
                if score > best_score:
                    best_score = score
                    best_match = (pathway_name, pathway_data)
            
            # Check description
            description = pathway_data.get('description', '').lower()
            if description and pathway_lower in description:
                score = len(pathway_lower) / len(description)
                if score > best_score:
                    best_score = score
                    best_match = (pathway_name, pathway_data)
        
        if best_match:
            pathway_name, pathway_data = best_match
            genes = pathway_data['genes']
            print(f"Found best match '{pathway_data.get('original_name', pathway_name)}': {len(genes)} genes")
            return genes
        
        print(f"No pathway found for: '{pathway_query}'")
        return []

    def get_functional_genes(self, function_name):
        """Get genes associated with a biological function"""
        matching_pathways = self.search_pathways_by_category(function_name)
        
        all_genes = []
        pathway_count = 0
        
        for pathway in matching_pathways:
            pathway_genes = self.get_pathway_genes(pathway['pathway_name'])
            if pathway_genes:
                all_genes.extend(pathway_genes)
                pathway_count += 1
        
        unique_genes = list(set(all_genes))
        print(f"Found {len(unique_genes)} unique genes from {pathway_count} pathways for function '{function_name}'")
        return unique_genes

    def list_available_categories(self):
        """List available functional categories"""
        return list(self.categories.keys())

    def get_library_stats(self):
        """Get statistics about loaded libraries"""
        library_stats = defaultdict(int)
        for pathway_data in self.pathways.values():
            library_stats[pathway_data['library']] += 1
        
        return dict(library_stats)

# Global database instance
pathway_db = None

def initialize_pathway_database(pathway_data_dir="pathway_data"):
    """Initialize the global pathway database"""
    global pathway_db
    if pathway_db is None:
        pathway_db = PathwayDatabase(pathway_data_dir)
    return pathway_db

def get_pathway_genes(pathway_name, action_type="pathway_filter"):
    """Get genes for a pathway"""
    global pathway_db
    if pathway_db is None:
        pathway_db = initialize_pathway_database()
    return pathway_db.get_pathway_genes(pathway_name, action_type)

def search_pathways_by_category(category):
    """Search pathways by category"""
    global pathway_db
    if pathway_db is None:
        pathway_db = initialize_pathway_database()
    return pathway_db.search_pathways_by_category(category)

def get_functional_genes(function_name):
    """Get genes for a biological function"""
    global pathway_db
    if pathway_db is None:
        pathway_db = initialize_pathway_database()
    return pathway_db.get_functional_genes(function_name)

def validate_pathway_command(action, target, value):
    """Validate pathway commands"""
    global pathway_db
    if pathway_db is None:
        pathway_db = initialize_pathway_database()
    
    if action == "pathway_filter":
        genes = pathway_db.get_pathway_genes(value)
        if not genes:
            return False, f"Pathway '{value}' not found. Try 'search {value} pathways' to see available options."
    
    elif action == "pathway_search":
        results = pathway_db.search_pathways_by_category(value)
        if not results:
            available_categories = pathway_db.list_available_categories()
            return False, f"No pathways found for '{value}'. Available categories: {', '.join(available_categories[:5])}, etc."
    
    elif action == "functional_filter":
        genes = pathway_db.get_functional_genes(value)
        if not genes:
            return False, f"No genes found for function '{value}'. Try searching for related pathways first."
    
    return True, ""

# Example usage
if __name__ == "__main__":
    db = PathwayDatabase("pathway_data")
    
    print("\n🧪 Testing pathway searches:")
    
    # Test cancer search
    cancer_results = db.search_pathways_by_category("cancer")
    print(f"Cancer pathways: {len(cancer_results)}")
    for i, result in enumerate(cancer_results[:3]):
        print(f"  {i+1}. {result['original_name'][:60]}... ({result['gene_count']} genes, {result['library']})")
    
    # Test immune search
    immune_results = db.search_pathways_by_category("immune")
    print(f"\nImmune pathways: {len(immune_results)}")
    for i, result in enumerate(immune_results[:3]):
        print(f"  {i+1}. {result['original_name'][:60]}... ({result['gene_count']} genes, {result['library']})")
    
    # Show library statistics
    print(f"\nLibrary statistics:")
    for library, count in db.get_library_stats().items():
        print(f"  {library}: {count} pathways")