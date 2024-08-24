from langchain_openai import ChatOpenAI
import os
from dotenv import load_dotenv
from openai import OpenAI
import nest_asyncio
nest_asyncio.apply()
import asyncio  # Add this import
import os
from dotenv import load_dotenv
import json
from APIs.combinedapi import PubMedProcessor
import datetime
from IPython.display import display, HTML
from difflib import unified_diff
import weave 
import weave
from weave import Evaluation
import json
import aiofiles
from difflib import unified_diff
import html







openai_client = OpenAI(api_key="OPENAI_API_KEY")
load_dotenv()
llm = ChatOpenAI()
llm = ChatOpenAI(api_key="OPENAI_API_KEY" )


chat_prompt = "You are ChatGPT, a large language model trained by OpenAI, based on the GPT-4 architecture and you cite the papers you use in your answers using Harvard Style."

working_hypothesis_prompt = """ 
# Scientific Rationale for 5-alpha reductase  in Huntington's Disease


## Target Information 
### Develop a scientific rationale for the following:
                             
    **Given target:**  5-alpa reductase
    **Given disease:** Huntington's disease
    **Given mode of action:** Inhibition of 5-alpha reductase can promote neuronal survival and reduce mutant huntingtin protein aggregation

## Task 1: Develop Scientific Rationale

### Working Hypothesis
- Detailed description of the idea
- Unmet medical need
- Suitability for combination therapy
- Predictive biomarkers
- Clinical relevance of existing biomarkers


"""

clinical_target_prompt = """ #Clinical Target Rationale for 5-alpha reductase in Huntington's Disease


## Target Information 
### Develop a scientific rationale for the following:
                           
    **Given target:** 5-alpa reductase
    **Given disease:** Huntington's disease
    **Given mode of action:** Inhibition of 5-alpha reductase can promote neuronal survival and reduce mutant huntingtin protein aggregation


 ### Clinical target rationale:
    - How relevant is the target location to the disease biology?
    - How is the target expression altered in human disease?
    - How is the target involved in the physiological process relevant to the disease?
    - Which phenotypes and genotypes were identified for the target?
    - How is the genetic link between the target and the disease?
    - Describe the evidence provided in clinics or by tools acting on the pathway where the target is involved.
    - Which kind of target modulation is required to treat the disease? """

Challenges_prompt_1 = """ 
#Challenges for the drug discovery program related to 5-alpha reductase  in Huntington's disease
                           
## Target Information 
### Develop a scientific rationale for the following:
                           
    **Given target:**  5-alpa reductase
    **Given disease:** Huntington's disease
    **Given mode of action:** Inhibition of 5-alpha reductase can promote neuronal survival and reduce mutant huntingtin protein aggregation


### Challenges:
- Check the following idea for details on small molecule compounds: Developing small molecule modulators or inhibitors of gamma secretase for Alzheimer's disease treatment.
- Is a 'information driven approach' (IDA) strategy based on available small molecules possible?
- Which small molecular modulators of the target known?
- Which inhibitors, antagonists, agonists, negative allosteric modulators (NAM), positive allosteric modulators (PAM) are required for target modulation in the given disease? 

"""
Challenges_prompt_2 = """

#Challenges for the drug discovery program related to 5-alpha reductase  in Huntington's disease
                           
## Target Information 
### Develop a scientific rationale for the following:
                           
    **Given target:**  5-alpa reductase
    **Given disease:** Huntington's disease
    **Given mode of action:** Inhibition of 5-alpha reductase can promote neuronal survival and reduce mutant huntingtin protein aggregation

- Which patients would respond the therapy?
- Is the proposed mode of action on the target desirable and commercially viable in a clinical setting?
- What are advantages and disadvantages of different therapeutic modalities (antibodies, small molecules, antisense oligonucleotides, PROTACs, molecular glue, peptide macrocycles, and so on) for tackling the target?

- Alternative indications:
- Describe alternative indication for modulators of the target and explain why.


"""
all_prompts = [working_hypothesis_prompt, clinical_target_prompt, Challenges_prompt_1, Challenges_prompt_2]

import json
from typing import Dict, Any
import openai
import asyncio
from typing import Dict, Any


async def pubmed_paperqa(query: str) -> Dict[str, Any]:
    """ Searches PubmedCentral for papers using a query
    and returns the most relevant chunks using paperQA"""

    max_attempts = 2
    max_results: int = 4
    pubmed_query = query
    doc_query = query
    email = "sanazkazemi@hotmail.com"
    print(f"pubmed_paperqa called with query: {query}, max_results: {max_results}")

    
    pubmed_instance = PubMedProcessor(email)

    # some error handeling in case the API call fails the algorithm will continue.
    results_dict = {}
    for attempt in range(max_attempts):
        try:
            results_dict = await pubmed_instance.full_process(pubmed_query, doc_query, max_results)
            break
        except Exception as e:
            print(f"Error in attempt {attempt+1}: {e}")

    if results_dict is None:
        print("All API call attempts failed. Continuing with empty results.")
        results_dict = {"error": "API calls failed", "results": []}


    # saving the refs to file because the LLM cant do it as it edits in paragraphs - has to be manually done

    file_path = 'RATT_refs.json'

    try:
        if os.path.exists(file_path) and os.path.getsize(file_path) > 0:
            with open(file_path, 'r') as f:
                existing_data = json.load(f)
        else:
            existing_data = []
    except json.JSONDecodeError:
        print("Error reading existing data. Starting with empty list.")

        existing_data = []

    existing_data.append(results_dict)

    with open(file_path, 'w') as f:
        json.dump(existing_data, f, indent=4)

                        
    return json.dumps(results_dict, indent=4)

# To run this in a Jupyter notebook cell:
# query = "Alzheimer's disease and gamma secretase"
# results = await pubmed_paperqa(query)
# print(results)


def generate_diff_html(text1, text2, fromfile='Original', tofile='Modified'):
    diff = unified_diff(text1.splitlines(keepends=True),
                        text2.splitlines(keepends=True),
                        fromfile=fromfile, tofile=tofile, n=3)
    
    html_output = ['''
    <style>
        .diff-container {
            font-family: monospace;
            white-space: pre-wrap;
            word-wrap: break-word;
            background-color: #f8f9fa;
            border: 1px solid #dee2e6;
            border-radius: 4px;
            padding: 10px;
            margin-bottom: 20px;
        }
        .diff-header {
            color: #6c757d;
            margin-bottom: 10px;
        }
        .diff-add {
            background-color: #e6ffec;
            color: #24292e;
        }
        .diff-sub {
            background-color: #ffebe9;
            color: #24292e;
        }
        .diff-line {
            display: block;
            margin-bottom: 0;
            padding: 2px 0;
        }
        .collapse-button {
            background-color: #007bff;
            color: white;
            border: none;
            padding: 5px 10px;
            margin-bottom: 10px;
            cursor: pointer;
            border-radius: 4px;
        }
        .hidden {
            display: none;
        }
    </style>
    <div class="diff-container">
    <button class="collapse-button" onclick="toggleDiff(this)">Collapse/Expand Diff</button>
    <div class="diff-content">
    ''']
    
    for line in diff:
        if line.startswith('---') or line.startswith('+++'):
            html_output.append(f'<div class="diff-header">{html.escape(line)}</div>')
        elif line.startswith('+'):
            html_output.append(f'<span class="diff-line diff-add">{html.escape(line)}</span>')
        elif line.startswith('-'):
            html_output.append(f'<span class="diff-line diff-sub">{html.escape(line)}</span>')
        else:
            html_output.append(f'<span class="diff-line">{html.escape(line)}</span>')
    
    html_output.append('''
    </div>
    </div>
    <script>
    function toggleDiff(button) {
        var content = button.nextElementSibling;
        if (content.style.display === "none") {
            content.style.display = "block";
            button.textContent = "Collapse Diff";
        } else {
            content.style.display = "none";
            button.textContent = "Expand Diff";
        }
    }
    </script>
    ''')
    
    return ''.join(html_output)




def pick_best_query(queries: list, question: str, answer: str) -> str:
    """picks the best query from the list of queries"""

    gpt_prompt ="""
    for the given question and answer, pick the best query 
    from the list of queries that you think is most relevant especially to the last few sentences of the answer.
    ## IMPORTANT:
    Just return the best query. Do not add any additional information.
    """

    best_query = openai.chat.completions.create(

     model = "gpt-4o-mini",
        messages=[
            {
                "role": "system",
                "content": "you are a scientific researcher, you are tasked with finding the best query to search for scientific papers on PubMed."            
        },
            {
                "role": "user",
                "content": f"##Instruction:{gpt_prompt}\n\n###Question: {question}\n\n###Question:{answer}\n\n##Queries: {[queries]}\n\n"
            }
        ],
        temperature=1 # here you can adjust the temperature to get more or less creative search terms
    ).choices[0].message.content

    return best_query




def get_query(question, answer, num_queries) -> str:

    """Generates queries to search for in Pubmed based on the question"""    
    query_prompt = """ You are a scientific researcher, 
                        you are tasked with finding the best query to search 
                        for scientific papers on PubMed.
                                                    
                            I want to verify the content correctness of the given answer especially the last few sentences.
                            Please summarize the content with the corresponding question.
                            This summarization will be used as a query to search with Bing search engine.
                            The query should be short but needs to be specific to promise Bing can find related knowledge or pages.
                            You can also use search syntax to make the query short and clear enough for the search engine to find relevant language data.
                            Try to make the query as relevant as possible to the last few setences of the the answer provided.
                            **IMPORTANT**
                            Just output the query directly. DO NOT add additional explanations or introducement in the answer unless you are asked to.

                        The following worked very well for me in the past in terms of generating the highest number of results use it as a guide:
                        ###Example:
                        "{Target}" AND "{Disease}" AND ("{relevant_keyword}" OR "{relevant_keyword_1} Or "{relevant_keyword_n}")" and so on.
                        ##IMPORTANT:
                        Just provide the query. Do not add any additional information.
                        DO NOT copy the given example"""
    

    queries = []

    for i in range(num_queries):
        try:
            query = openai.chat.completions.create(
                model="gpt-4o-mini",
                messages=[
                    {
                        "role": "system",
                        "content": "You are ChatGPT, a large language model trained by OpenAI, based on the GPT-4 architecture."

                    },
                    {
                        "role": "user",
                        "content": f"##Question: {question}\n\n##Content: {answer}\n\n##Instruction: {query_prompt}"
                    }
                ],
                temperature=1
            ).choices[0].message.content

            print(f"query {i}: {query}")

        except Exception as e:
            print(f"error {e}")
        queries.append(query)

    best_query = pick_best_query(queries,answer, question)
    print(f"best query: {best_query}")

    return best_query




async def main(question: str, answer, num_queries: int):
    """Main function to get the best query for the question"""
    
    
    search_query = get_query(question, answer, num_queries)
    
    # Remove only the outermost single quotes if they exist otherwise doesnt work - not elegant but works
    if search_query.startswith("'") and search_query.endswith("'"):
        cleaned_query = search_query[1:-1]
    else:
        cleaned_query = search_query
    
    # Replace escaped single quotes with regular single quotes
    cleaned_query = cleaned_query.replace("\\'", "'")
    
    print(f"Cleaned search query: {cleaned_query}")
    
    results = await pubmed_paperqa(cleaned_query)

    print(results)  # This is the final output


    return results

num_agents = 1 
num_steps = 3
final_output_mode = 'final_step_only'


def COT_agent(question):
    """Generates a chain of thought answer for comparison to RATT.
    question: str: the prompt to answer
    draft_prompt: str: the prompt to generate the draft
    system_prompt: str: the prompt to generate the system message"""

    draft_prompt = '''
IMPORTANT:
Try to answer this question/instruction with step-by-step thoughts and make the answer more structured.
Use `\n\n` to split the answer into several paragraphs.
Just respond to the instruction directly. DO NOT add additional explanations or introducement in the answer unless you are asked to.
'''

    # Loop to generate different initial answers
    COT_draft = openai.chat.completions.create(
         model="gpt-4o-mini",
            messages=[
                {
                    "role": "system",
                    "content": chat_prompt # should this be scientific rationale prompt or something less specific?
                },
                {
                    "role": "user",
                    "content": question + draft_prompt
                }
            ],
            temperature=0.5
        ).choices[0].message.content

    return COT_draft

from datetime import datetime


def split_draft(draft, split_char='\n\n'):
    # split_char: '\n\n'
    draft_paragraphs = draft.split(split_char)
    # print(f"The draft answer has {len(draft_paragraphs)}")
    return draft_paragraphs


def get_revise_answer(question, answer, retrieved_data):
    revise_prompt = '''
I want to revise the answer according to retrieved related text of the question from Pubmed central. You will receive the summary, full text, relevance score to the query, and the full reference. If you use any of the given articles, you **must** reference using the Harvard style with DOI added and add the full citation to the top of the answer. 
**DO NOT REMOVE OR ALTER ANY CITATIONS FROM THE TOP OF THE ANSWER.**
You need to check whether the answer is correct.
If you find some errors in the answer, revise the answer to make it better.
If you find some necessary details are ignored, add it to make the answer more plausible according to the related text.
If you find that a part of the answer is correct and does not require any additional details, maintain that part of the answer unchanged. Directly output the original content of that part without any modifications.
**IMPORTANT**
Try to keep the structure (multiple paragraphs with its subtitles) in the revised answer and make it more structual for understanding.
Split the paragraphs with `\n\n` characters.
Just output the revised answer directly. DO NOT add additional explanations or annoucement in the revised answer unless you are asked to.
'''
    revised_answer = openai.chat.completions.create(
        model="gpt-4o-mini",
        messages=[
            {
                "role": "system",
                "content": chat_prompt
            },
            {
                "role": "user",
                "content": f"##Pubmed central retrieved articles: {retrieved_data}\n\n##Question: {question}\n\n##previous Answer: {answer}\n\n##Instruction: {revise_prompt}"
            }
        ],
        temperature=0.5
    ).choices[0].message.content
    
    
    return revised_answer



async def RAG(question, draft_paragraphs):
    """ args:
    question: str: the prompt to answer
    draft_paragraphs: list: the list of paragraphs from the initial n drafts
    """
    answer = ""

    for i, paragraph in enumerate(draft_paragraphs):
        answer += '\n\n' + paragraph

        api_response = await main(question, answer, num_queries=2)  # Now using the entire answer instead of just the paragraph

        revised_answer = get_revise_answer(question, answer, api_response)  # Using the entire answer
        if revised_answer != answer:
            diff_html = generate_diff_html(answer, revised_answer)
            display(HTML(diff_html))
            answer = revised_answer
        
        print(f"Completed iteration {i+1}/{len(draft_paragraphs)}")

        print('+'* 80 + '\n\n')
        print(f"RESULT OF PUBMED API:\n{answer}")
    
    return answer



async def get_draft_tot_initial(question: str, num_agents: int):
    """Generates initial answers from multiple agents for comparison"""
    draft_prompt = """
            IMPORTANT:
            Try to answer this question/instruction with step-by-step thoughts and make the answer more structured.
            Use `\n\n` to split the answer into several paragraphs.
            Just respond to the instruction directly. DO NOT add additional explanations or introducement in the answer unless you are asked to.
            """

    refine_prompt = """
            Maintaining *ALL* citations and references and referencing the answers provided by all agents, synthesize a more detailed and comprehensive response by integrating all relevant details from these answers. 
            Ensure logical coherence and provide ONLY THE MERGED ANSWER AND CITATIONS AS THE OUTPUT, omitting any discussion of the comparison process or analytical thoughts."""

    agent_drafts = []
    for i in range(num_agents):
        draft = openai.chat.completions.create(
            model="gpt-4o-mini",
            messages=[
                {
                    "role": "system",
                    "content": chat_prompt
                },
                {
                    "role": "user",
                    "content": question + draft_prompt
                }
            ],
            temperature=0.5
        ).choices[0].message.content
        print(f"####################draft {i}: {draft}########################################")

        print("Processing draft...")
        draft_paragraphs = split_draft(draft)

        draft_modified = await RAG(question, draft_paragraphs)

        agent_drafts.append(f"Agent{i+1}: {draft_modified}")

        print(f"[INFO] Agent{i + 1}/{num_agents} retrieved draft...")

        agent_input = '\n\n'.join(agent_drafts) + '\n\n' + refine_prompt

        final_draft = openai.chat.completions.create(
            model="gpt-4o-mini",
            messages=[
                {
                    "role": "system",
                    "content": chat_prompt
                },
                {
                    "role": "user",
                    "content": agent_input
                }
            ],
            temperature=0.5
        ).choices[0].message.content

        print(f"{datetime.now()} - Final draft: {final_draft}")


    return final_draft

# FIX BELOW FUNCTION NOT USING PREVIOUS ANSWER ?? ALSO STRING CHAR LIST WITHIN LIST - SHOULD BE STRING


async def get_draft_tot(question, previous_answer, num_agents):

    draft_prompt = f""" Base your response on the provided question and the previous answer. Expand the answer by adding more details to enhance its comprehensiveness. Ensure that the expansion maintains logical coherence and enriches the details, making the response more thorough and well-structured.
        Question: {question}
        Previous Answer: {previous_answer}
        IMPORTANT:
        DO NOT REMOVE ANY CITATIONS OR REFERENCES IN THE ANSWER if you use any citations or references in the answer - you must reference in Harvard style + the DOI in the answer and add the citation to the top of the answer.
        Answer the full question with step-by-step thoughts and make the answer more structural.
        Use `\n\n` to split the answer into several paragraphs.
        Just respond to the instruction directly. DO NOT add additional explanations or introducement in the answer unless you are asked to. """

    refine_prompt = """Maintaining *ALL* citations and references and referencing the answers provided by all agents, synthesize a more detailed and comprehensive response by integrating all relevant details from these answers. 
            Ensure logical coherence and provide ONLY THE MERGED ANSWER AND CITATIONS AS THE OUTPUT, omitting any discussion of the comparison process or analytical thoughts."""

    agents_drafts = []
    for i in range(num_agents):
        draft = openai.chat.completions.create(
            model="gpt-4o-mini",
            messages=[
                {
                    "role": "system",
                    "content": chat_prompt
                },
                {
                    "role": "user",
                    "content": draft_prompt
                }
            ],
            temperature=0.5
        ).choices[0].message.content

        draft_paragraphs = split_draft(draft)

        draft_modified = await RAG(question, draft_paragraphs)

        agents_drafts.append(f"Agent{i+1}: {draft_modified}")
    
    agents_input = '\n\n'.join(agents_drafts) + '\n\n' + refine_prompt

    final_draft_raw = openai.chat.completions.create(
        model="gpt-4o-mini",
        messages=[
            {
                "role": "system",
                "content": chat_prompt
            },
            {
                "role": "user",
                "content": agents_input
            }
        ],
        temperature=0.5
    ).choices[0].message.content

    print(f"##########Final draft raw #########################: {final_draft_raw}...")

    revise_prompt = """
            Based on the original answer and an additional supplementary answer, generate a response that is richer in detail and logically coherent. MAINTAIN ALL REFERENCES AND CITATIONS IN THE ANSWER. Do not remove ANY citations.
            Review the original answer:
        1. If any part of the answer is correct and requires no further details, retain that portion unchanged and output it directly as it is.
        2. For parts that may be improved or lack necessary details, enhance them by integrating information from the supplementary answer to make the response more comprehensive and accurate, reference and cite the full reference at the top where neccessary.
        3. If you identify any errors within the answers, correct these errors while ensuring that the revised content remains logically coherent - check this against the retreived articles.
        Original Answer: {previous_answer}
        Supplementary Answer: {final_draft_raw}

        **IMPORTANT**
        Ensure the revised answer maintains a structured format with paragraphs and subtitles for better clarity. 
        Separate the paragraphs with `\n\n` characters. Output only the enhanced answer directly and the references and citations if any, without any extra explanations or announcements unless specifically requested."""
    
    final_draft = openai.chat.completions.create(
        model="gpt-4o-mini",
        messages=[
            {
                "role": "system",
                "content": chat_prompt
            },
            {
                "role": "user",
                "content": revise_prompt
            }
        ],
        temperature=0.5
    ).choices[0].message.content

    return final_draft




async def ratt(question, num_agents):
    step_num = num_steps
    print(f"{datetime.now()} [INFO] Retrieving Step 1 draft...")

    draft = await get_draft_tot_initial(question,num_agents)
    
    print(f"{datetime.now()} [INFO] Step 1 draft returned")
    print(f"##################### DRAFT #######################")
    print(draft)
    print(f"#####################  END  #######################")

    print(f"{datetime.now()} [INFO] Processing draft...")
    draft_paragraphs = split_draft(draft)
    print(f"{datetime.now()} [INFO] Draft split into {len(draft_paragraphs)} parts")

    answer_first_state = await RAG(question, draft_paragraphs)

    previous_answer = answer_first_state

    each_step_drafts = [f"Step 1 \n: {previous_answer}"]

    for iteration in range(1, step_num):
        print(f"{datetime.now()} [INFO] Retrieving Step {iteration + 1} draft...")
        draft = await get_draft_tot(question, previous_answer, num_agents=num_agents)
        print(f"{datetime.now()} [INFO] Step {iteration + 1} draft returned")
        print(f"##################### DRAFT #######################")
        print(draft)
        print(f"#####################  END  #######################")

        print(f"{datetime.now()} [INFO] Processing draft...")
        draft_paragraphs = split_draft(draft)
        print(f"{datetime.now()} [INFO] Draft split into {len(draft_paragraphs)} parts")

        # filtered_paragraphs = filter_paragraphs(draft_paragraphs, iteration, step_num)
        final_answer = await RAG(question, draft_paragraphs)

        each_step_drafts.append(f"Step {iteration + 1} \n: {final_answer}")

        # Update previous_answer for the current iteration's response
        previous_answer = final_answer

    draft_cot = COT_agent(question) # for comparison

    if final_output_mode == 'combine_each_step':
        final_draft = '\n\n'.join(each_step_drafts)
        refine_prompt = f"""
            Maintaining *ALL* citations and references and referencing the answers provided by all agents, synthesize a more detailed and comprehensive response by integrating all relevant details from these answers. 
            Ensure logical coherence and provide ONLY THE MERGED ANSWER AND CITATIONS AS THE OUTPUT, omitting any discussion of the comparison process or analytical thoughts."""
        previous_answer = openai.chat.completions.create(
            model="gpt-4o-mini",
            messages=[
                {
                    "role": "system",
                    "content": chat_prompt
                },
                {
                    "role": "user",
                    "content": final_draft + '\n\n' + refine_prompt
                }
            ],
            temperature=0.5 ## ? why
        ).choices[0].message.content

    return draft_cot, previous_answer




@weave.op()
async def main_ratt(question, num_agents, output_file="text_ratt/EGFR/EGFR.json"):
    answer_cot, answer_ratt = await ratt(question, num_agents)

    result = {
        "prompt": question,
        "output_cot": answer_cot, 
        "output_ratt": answer_ratt
    }

    # Write the result to a file
    async with aiofiles.open(output_file, mode='a') as f:
        await f.write(json.dumps(result) + "\n")

    return result

    
weave.init('output_ratt_EGFR')
result = asyncio.run(main_ratt("What is the role of EGFR in cancer?", num_agents=1))

# weave.init('output_ratt_bace')
# for prompt in all_prompts:
#     print(f"**************************{prompt}****************************")
#     results = asyncio.run(main_ratt(prompt, num_agents=1))

