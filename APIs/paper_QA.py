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


    prompts = PromptCollection(
            summary_prompt = (
            "Extract key information from the excerpt below that could be useful for other AI agents.\n\n"
            "Excerpt from {citation}\n\n----\n\n{text}\n\n----\n\n"
            "Focus: {question}\n\n"
            "Do not form a summary or answer the question directly."
            "Include specific numbers, equations, or direct quotes (marked with quotation marks) when relevant. "
            'Reply "Not applicable" if the excerpt contains no useful information. '
            "At the end, provide an integer score from 1-10 on a newline indicating usefulness of the information. "
            "Do not explain your score."
            "\n\nKey Information Summary ({summary_length}):"
        ),

        summary_json_prompt = (
            "Excerpt from {citation}\n\n----\n\n{text}\n\n----\n\nFocus: {question}\n\n"
        ),

        qa_prompt = (
            "Compile relevant information from the context below.\n\n"
            "Context (with usefulness scores):\n\n{context}\n\n----\n\n"
            "Focus: {question}\n\n"
            "Provide a comprehensive compilation of relevant information from the context. "
            "Do not form conclusions or direct answers. "
            "For each piece of information, indicate the source via citation keys at the end of sentences, "
            "like (Author, Year, p. 3-4). Only cite from the provided context. "
            "Organize the information in a clear, structured manner. "
            "If there are conflicting data or perspectives, present them objectively. "
            "Include relevant quotes if present.\n\n"
            "Information Compilation ({answer_length}):"
        ),

        select_paper_prompt = (
            "Select papers that may help answer the question below. "
            "Papers are listed as $KEY: $PAPER_INFO. "
            "Return a list of keys, separated by commas. "
            'Return "None", if no papers are applicable. '
            "Choose papers that are relevant, from reputable sources, and timely "
            "(if the question requires timely information).\n\n"
            "Question: {question}\n\n"
            "Papers: {papers}\n\n"
            "Selected keys:"
        ),

        citation_prompt = (
            "Provide the citation for the following text in Harvard citation format. "
            "Do not write an introductory sentence. "
            "If reporting date accessed, the current year is 2024\n\n"
            "{text}\n\n"
            "Citation:"
        ),

        default_system_prompt = (
            "Compile information in a direct and specific manner. "
            "Your output will be used by other AI agents, so focus on facts and data. "
            "Define any ambiguous terms or acronyms. "
            "Do not form conclusions or provide direct answers to questions."
        ),

        summary_json_system_prompt = """\
        Provide a summary of the key information from the excerpt that could be useful for other AI agents. Respond with the following JSON format:
        {
        "key_information": "...",
        "usefulness_score": "..."
        }
        where `key_information` is relevant data from text - {summary_length} words and `usefulness_score` is the potential usefulness of this information for other AI agents (out of 10).
        """
    )

    docs = Docs(prompts=prompts)
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

docs = add_docs("/Users/sanazkazeminia/Documents/LLM_Agent/Target-Assessment-Agent-/pmc_papers/")
answer = query_docs("Machine Learning AND COVID-19", docs)
context_for_qa = answer.context
print(context_for_qa)


