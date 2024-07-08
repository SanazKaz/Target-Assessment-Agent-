import xml.etree.ElementTree as ET
import json
from collections import defaultdict

def generate_summary(drug_data):
    summary = f"Drug {drug_data['name']} (DrugBank ID: {drug_data['drugbank-id']}) "
    if drug_data.get('indication'):
        summary += f"is indicated for {drug_data['indication']} "
    if drug_data.get('mechanism-of-action'):
        summary += f"It acts by {drug_data['mechanism-of-action']}"
    return summary.strip()

def extract_drug_info(xml_file, output_file):
    target_elements = ['drugbank-id', 'name', 'description', 'indication', 'mechanism-of-action',
                       'pharmacodynamics', 'toxicity', 'targets', 'enzymes', 'carriers', 'transporters']

    drugs = defaultdict(lambda: defaultdict(str))
    ns = {'db': 'http://www.drugbank.ca'}

    for event, elem in ET.iterparse(xml_file, events=('end',)):
        if elem.tag == '{http://www.drugbank.ca}drug':
            drug_id = elem.find('db:drugbank-id', ns).text
            for el in target_elements:
                el_with_ns = f"{{http://www.drugbank.ca}}{el}"
                if elem.find(el_with_ns, ns) is not None:
                    if el in ['targets', 'enzymes', 'carriers', 'transporters']:
                        items = elem.find(el_with_ns, ns).findall('.//db:name', ns)
                        drugs[drug_id][el] = [item.text for item in items if item.text]
                    else:
                        drugs[drug_id][el] = elem.find(el_with_ns, ns).text or ""
            
            # Generate summary without truncation
            drugs[drug_id]['summary'] = generate_summary(drugs[drug_id])
            
            elem.clear()

    # Convert defaultdict to list of drugs
    drug_list = [dict(drug_data) for drug_data in drugs.values()]

    with open(output_file, 'w', encoding='utf-8') as jsonfile:
        json.dump(drug_list, jsonfile, ensure_ascii=False, indent=2)

# Usage
extract_drug_info('/Users/sanazkazeminia/Documents/LLM_Agent/drugbank_data/full database.xml', '/Users/sanazkazeminia/Documents/LLM_Agent/drugbank_data/drug_target_info.json')
print("Improved drug information has been extracted to improved_drug_info.json")
