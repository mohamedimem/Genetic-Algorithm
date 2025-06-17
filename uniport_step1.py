from unipressed import UniprotkbClient
from typing import Set, Dict, List, Any

class UniProtGOExtractor:
    """
    Step 1: UniProt GO:CC Annotation Research and Extraction
    
    Purpose: Extract official Gene Ontology Cellular Component (GO:CC) 
    annotations for a given human gene from UniProt using the Unipressed library.
    """
    
    def __init__(self):
        self.client = UniprotkbClient
        
    def step1_fetch_go_annotations(self, gene_name: str) -> Set[str]:
        """
        Extract GO:CC annotations for a human gene from UniProt.
        
        Logic Breakdown (step-by-step):
        1. Query the UniProt REST API with:
           * Gene name (e.g., "POLE3")
           * Filter for:
              * Organism: Human (organism_id=9606)
              * Reviewed entries (Swiss-Prot)
        2. Parse the JSON response and extract GO:CC annotations from:
           * goAnnotations field with aspect == "C"
           * uniProtKBCrossReferences where id starts with "GO:" and properties.aspect == "C"
        3. Collect and return all matching GO:CC terms as a set (to avoid duplicates).
        
        Args:
            gene_name (str): Human gene name (e.g., "POLE3")
            
        Returns:
            Set[str]: Set of GO:CC terms. Example: {'GO:0005634', 'GO:0005737', 'GO:0005739'}
        """
        
        print(f"🔍 Searching UniProt for gene: {gene_name}")
        
        try:
            # Step 1: Search UniProt using Unipressed with proper query structure
            search_results = self.client.search(
                query={
                    "and_": [
                        {"gene": gene_name},
                        {"organism_id": "9606"},  # Human
                        {"reviewed": "true"}      # Swiss-Prot only
                    ]
                }
            )
            
            # Step 2: Process each record found
            go_cc_terms = set()
            record_count = 0
            
            for record in search_results.each_record():
                record_count += 1
                print(f"📋 Processing record {record_count}: {record.get('primaryAccession', 'Unknown')}")
                
                # Step 3: Extract GO:CC annotations using both methods
                # Method 1: Extract from goAnnotations field
                go_cc_terms.update(self._extract_from_go_annotations(record))
                
                # Method 2: Extract from uniProtKBCrossReferences field  
                go_cc_terms.update(self._extract_from_cross_references(record))
            
            
            return go_cc_terms
            
        except Exception as e:
            print(f"❌ Error during UniProt search: {e}")
            return set()
    
    def _extract_from_go_annotations(self, record: Dict[str, Any]) -> Set[str]:
        """
        Extract GO:CC terms from the goAnnotations field.
        Looking for annotations where aspect == "C"
        """
        go_cc_terms = set()
        
        go_annotations = record.get('goAnnotations', [])
        if go_annotations:
            print(f"  🔍 Processing {len(go_annotations)} GO annotations...")
        
        for annotation in go_annotations:
            aspect = annotation.get('aspect')
            go_id = annotation.get('goId') 
            
            if aspect == 'C':  # Cellular Component
                if go_id and go_id.startswith('GO:'):
                    go_cc_terms.add(go_id)
                    print(f"  ✅ Found GO:CC from goAnnotations: {go_id}")
        
        return go_cc_terms
    
    def _extract_from_cross_references(self, record: Dict[str, Any]) -> Set[str]:
        """
        Extract GO:CC terms from uniProtKBCrossReferences field.
        Looking for entries where database == "GO" and properties indicate aspect == "C"
        """
        go_cc_terms = set()
        
        cross_refs = record.get('uniProtKBCrossReferences', [])
        if cross_refs:
            print(f"  🔍 Processing {len(cross_refs)} cross-references...")
        
        go_refs_count = 0
        for cross_ref in cross_refs:
            if cross_ref.get('database') == 'GO':
                go_refs_count += 1
                go_id = cross_ref.get('id')
                properties = cross_ref.get('properties', [])
                
                if go_id and go_id.startswith('GO:'):
                    # Check properties for aspect information
                    is_cellular_component = False
                    
                    for prop in properties:
                        if prop.get('key') == 'GoTerm':
                            # Check if this is a cellular component based on term description
                            term_desc = prop.get('value', '').lower()
                            if any(cc_indicator in term_desc for cc_indicator in 
                                   ['cellular', 'component', 'compartment', 'organelle', 
                                    'membrane', 'nucleus', 'cytoplasm', 'mitochondrial',
                                    'ribosome', 'endoplasmic', 'golgi', 'peroxisome',
                                    'vacuole', 'lysosome', 'centrosome']):
                                is_cellular_component = True
                                break
                        elif prop.get('key') == 'GoEvidenceType':
                            # Some implementations may have aspect info here
                            evidence = prop.get('value', '')
                            if 'C' in evidence or 'cellular' in evidence.lower():
                                is_cellular_component = True
                                break
                    
                    if is_cellular_component:
                        go_cc_terms.add(go_id)
                        print(f"  ✅ Found GO:CC from crossReferences: {go_id}")
        
        if go_refs_count > 0:
            print(f"  📊 Found {go_refs_count} GO cross-references total")
        
        return go_cc_terms
    
    def get_detailed_annotations(self, gene_name: str) -> Dict[str, Any]:
        """
        Optional method: Get detailed GO annotations with metadata
        Returns full annotation details for debugging/analysis purposes
        """
        print(f"🔍 Getting detailed annotations for gene: {gene_name}")
        
        try:
            search_results = self.client.search(
                query={
                    "and_": [
                        {"gene": gene_name},
                        {"organism_id": "9606"},
                        {"reviewed": "true"}
                    ]
                }
            )
            
            detailed_results = {
                'gene_name': gene_name,
                'records': [],
                'go_cc_annotations': [],
                'go_cc_cross_refs': []
            }
            
            for record in search_results.each_record():
                record_info = {
                    'accession': record.get('primaryAccession'),
                    'protein_name': record.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', 'Unknown'),
                    'gene_names': [gene.get('geneName', {}).get('value') for gene in record.get('genes', [])],
                }
                detailed_results['records'].append(record_info)
                
                # Collect detailed GO annotations
                for annotation in record.get('goAnnotations', []):
                    if annotation.get('aspect') == 'C':
                        detailed_results['go_cc_annotations'].append(annotation)
                
                # Collect GO cross-references
                for cross_ref in record.get('uniProtKBCrossReferences', []):
                    if cross_ref.get('database') == 'GO':
                        detailed_results['go_cc_cross_refs'].append(cross_ref)
            
            return detailed_results
            
        except Exception as e:
            print(f"❌ Error getting detailed annotations: {e}")
            return {}
    
    def validate_go_terms(self, go_terms: Set[str]) -> Dict[str, bool]:
        """
        Optional method: Validate that extracted GO terms are properly formatted
        Returns a dictionary mapping each term to its validity status
        """
        validation_results = {}
        
        for term in go_terms:
            # GO terms should match pattern: GO:XXXXXXX (where X is a digit)
            is_valid = (term.startswith('GO:') and 
                       len(term) == 10 and 
                       term[3:].isdigit())
            validation_results[term] = is_valid
            
            if not is_valid:
                print(f"⚠️ Invalid GO term format: {term}")
        
        return validation_results

# Installation Instructions:
# pip install unipressed

# Example usage and testing
if __name__ == "__main__":
    # Initialize the extractor
    extractor = UniProtGOExtractor()
 
    go_cc_terms = extractor.step1_fetch_go_annotations("POLE3")
    
    print("\n📋 Results Summary:")
    print(f"Gene: POLE3")
    print(f"GO:CC Terms Found: {len(go_cc_terms)}")
    

    """ 💾 Return type: <class 'set'>
    Example return: {'GO:0006334', 'GO:0005721', 'GO:0008623', 'GO:0006261', 'GO:0000122', 'GO:0140672', 'GO:0006338', 'GO:0031490', 'GO:0006974', 
    'GO:0031507', 'GO:0006272', 'GO:0005634', 'GO:0006275'} 
    print(f"\n💾 Return type: {type(go_cc_terms)}")
        print(f"Example return: {go_cc_terms}")
    """
   
    print("\n🎉 Step 1 Implementation Complete!")
    print("Ready for Step 2: GO Term Enrichment Analysis")