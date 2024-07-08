import ollama
import json
import re

# Read the JSON file
with open('/Users/sanazkazeminia/Documents/LLM_Agent/drugbank_data/drug_target_info.json', 'r', encoding='utf-8') as jsonfile:
    drugs = json.load(jsonfile)

if drugs:
    print(f"Loaded {len(drugs)} drugs")
else:
    print("File not loaded")

def extract_drug_info(drug_name):
    """
    Extract information for a specific drug from the loaded data.
    """
    drug_name_lower = drug_name.lower()
    for drug in drugs:
        if drug['name'].lower() == drug_name_lower:
            return drug
    return None

def extract_drug_names_from_prompt(prompt):
    """
    Extract potential drug names from the given prompt.
    """
    # This regex will match words that start with a letter and are at least 3 characters long
    words = re.findall(r'\b[a-zA-Z][a-zA-Z\-]{2,}\b', prompt)
    # Filter out common words that are unlikely to be drug names
    common_words = set(['the', 'and', 'for', 'with', 'what', 'how', 'why', 'does', 'is', 'are', 'was', 'were'])
    return [word for word in words if word.lower() not in common_words]

def enrich_prompt_with_drug_info(prompt, drug_info):
    """
    Enrich the original prompt with extracted drug information.
    """
    if drug_info:
        additional_info = f"\n\nAdditional information about {drug_info['name']}:\n"
        additional_info += f"Description: {drug_info['description']}\n"
        additional_info += f"Targets: {', '.join(drug_info['targets'])}\n"
        additional_info += f"Mechanism of Action: {drug_info['mechanism-of-action']}\n"
        return prompt + additional_info
    return prompt

def get_llm_response(prompt):
    """
    Get a response from the LLM using the provided prompt.
    """
    response = ollama.chat(model='llama3', messages=[
        {
            'role': 'user',
            'content': prompt,
        },
    ])
    return response['message']['content']

def main():
    user_prompt = input("Enter your prompt: ")
    print("\n<agent> Original Prompt:")
    print(user_prompt)
    
    potential_drug_names = extract_drug_names_from_prompt(user_prompt)
    print(f"\n<agent> Potential drug names extracted: {potential_drug_names}")
    
    enriched_prompt = user_prompt
    found_drugs = []

    for drug_name in potential_drug_names:
        drug_info = extract_drug_info(drug_name)
        if drug_info:
            enriched_prompt = enrich_prompt_with_drug_info(enriched_prompt, drug_info)
            found_drugs.append(drug_name)
        else:
            print(f"<agent> No information found for {drug_name}")

    if found_drugs:
        print(f"\n<agent> Found information for the following drugs: {', '.join(found_drugs)}")
        print("\n<agent> Enriched Prompt:")
        print(enriched_prompt)
    else:
        print("\n<agent> No drug information found to enrich the prompt.")

    print("\n<agent> Generating a response...")
    response = get_llm_response(enriched_prompt)
    print("\n<agent> Response:")
    print(response)

if __name__ == "__main__":
    main()