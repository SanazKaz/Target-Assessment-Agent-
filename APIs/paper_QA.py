from paperqa import Docs, Answer, PromptCollection
from paperqa.prompts import select_paper_prompt, citation_prompt
import os
import json
import openai
from dotenv import load_dotenv


load_dotenv()
openai_api_key = os.getenv("OPENAI_API_KEY")


def add_docs(directory_path):
    """Adds documents to the document store
    and returns the document store object"""


    prompts = PromptCollection()

    docs = Docs(prompts=prompts, llm="gpt-4o-mini")
    print(select_paper_prompt)
    
    my_docs = []
    
    for root, dirs, files in os.walk(directory_path):
        for file in files:
            if file.endswith(".txt"):
                path = os.path.join(root, file)
                my_docs.append(path)
    
    print(f"these are your documents: {my_docs}")
    
    for d in my_docs:
        print(f"Adding {d} to the document store...")
        try:
            docs.add(d)
        except Exception as e:
            print(f"Error adding {d}: {e}")
    
    return docs

def query_docs(query, docs):
    """Queries the document store and returns the answer using paperqa"""
    answer = docs.query(query)
    
    return answer

# for text in docs.texts: # accesses ALL the chunks not just relevant ones
#     print(text.name, text.text)

docs = add_docs("/Users/sanazkazeminia/Documents/LLM_Agent/Target-Assessment-Agent-/pmc_papers/")
answer_1 = query_docs("relevance of dopamine 1 receptor to parkinsons disease", docs)
# print(f"answer: {answer_1.answer}")
context_answer = answer_1.context
# print(f"context: {context_answer}")
evidence = docs.get_evidence(answer_1, k=10, max_sources=10, detailed_citations=True)

# print(f"evidence: {evidence}")
# for text in docs.texts: # accesses ALL the chunks not just relevant ones
#     print(text.name, text.text)

# for context in evidence.contexts:
#     print(context.text.name, context.text.text)
    


dictionary_for_llm = {}

for context in evidence.contexts:
    # The key will be the context summary
    summary = context.context  # This is the summary/context we got from the LLM

    source_info = {
        "chunk_id": context.text.name,
        "full_citation": context.text.doc.citation,
        }



    # The value will be a dictionary containing the original text, source, and relevance score
    chunk_info = {
        "original_text": context.text.text,
        "source": source_info,
        "relevance_score": context.score  # Adding the relevance score
    }

    # Add to our main dictionary, using the summary as the key
    dictionary_for_llm[summary] = chunk_info

# print("keys", dictionary_for_llm.keys())
# print("values", dictionary_for_llm.values())
# print(len(dictionary_for_llm))

print("dictionary_for_llm", dictionary_for_llm)

with open('dictionary_for_llm.json', 'w') as f:
    json.dump(dictionary_for_llm, f, indent=4)