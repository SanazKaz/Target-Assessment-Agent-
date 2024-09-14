from IPython.display import display, HTML
from difflib import unified_diff
from langchain_openai import ChatOpenAI
from openai import OpenAI
import os
import openai
from dotenv import load_dotenv
import json
import datetime
import weave 
import weave
import json
from multiprocessing import Process, Queue
from langchain_community.document_loaders import PyPDFLoader
from langchain_openai import OpenAIEmbeddings
from langchain.schema import Document
from llama_index.core.node_parser import SentenceWindowNodeParser
from llama_index.core.postprocessor import MetadataReplacementPostProcessor
from llama_index.llms.openai import OpenAI
from llama_index.core import Settings
from langchain_community.document_loaders import PyPDFLoader
from llama_index.core import Document, VectorStoreIndex, Settings
from langchain_openai import OpenAIEmbeddings, OpenAI
from llama_index.core.node_parser import SentenceWindowNodeParser
from llama_index.core.postprocessor import MetadataReplacementPostProcessor
from llama_index.llms.openai import OpenAI
from llama_index.core import (
    load_index_from_storage,
    StorageContext
)
query_engine = None

num_agents = 1
num_steps = 3
final_output_mode = 'final_step_only'


# Directory containing your PDFs
pdf_directory = "/Users/sanazkazeminia/Documents/LLM_Agent/new_set_articles"

# for storing th pdfs
all_docs = []


# all prompts single and combined

chat_prompt = "You are a highly trained scientist and an expert in scientific writing you cite the papers you use in your answers using the name of the PDF provided to you."

working_hypothesis_prompt = """ 
# Scientific Rationale gamma secretase in Alzheimer's Disease


## Target Information 
### Develop a scientific rationale for the following:
                            
    **Given target:** gamma secretase
    **Given disease:** Alzheimer's disease
    **Given mode of action:** Inhibition of gamma secretase can reduce the production of amyloid-beta peptides and prevent the formation of amyloid plaques in the brain of Alzheimer's disease patients.
    
## Task 1: Develop Scientific Rationale

### Working Hypothesis
- Detailed description of the idea
- Unmet medical need
- Suitability for combination therapy
- Predictive biomarkers
- Clinical relevance of existing biomarkers


"""

clinical_target_prompt = """ #Clinical Target Rationale for gamma sec in Alzheimer's Disease

## Target Information 
### Develop a scientific rationale for the following:
                            
    **Given target:** gamma secretase
    **Given disease:** Alzheimer's disease
    **Given mode of action:** Inhibition of gamma secretase can reduce the production of amyloid-beta peptides and prevent the formation of amyloid plaques in the brain of Alzheimer's disease patients.
    
### Clinical target rationale:
    - How relevant is the target location to the disease biology?
    - How is the target expression altered in human disease?
    - How is the target involved in the physiological process relevant to the disease?
    - Which phenotypes and genotypes were identified for the target?
    - How is the genetic link between the target and the disease?
    - Describe the evidence provided in clinics or by tools acting on the pathway where the target is involved.
    - Which kind of target modulation is required to treat the disease? """


Challenges_prompt_1 = """ 
#Challenges for the drug discovery program related to gamma secretase in Alzheimer's disease
                        
## Target Information 
### Develop a scientific rationale for the following:
                            
    **Given target:** gamma secretase
    **Given disease:** Alzheimer's disease
    **Given mode of action:** Inhibition of gamma secretase can reduce the production of amyloid-beta peptides and prevent the formation of amyloid plaques in the brain of Alzheimer's disease patients.
    

### Challenges:
- Check the following idea for details on small molecule compounds: Developing small molecule modulators or inhibitors of gamma secretase for Alzheimer's disease treatment.
- Is a 'information driven approach' (IDA) strategy based on available small molecules possible?
- Which small molecular modulators of the target known?
- Which inhibitors, antagonists, agonists, negative allosteric modulators (NAM), positive allosteric modulators (PAM) are required for target modulation in the given disease? 

"""
Challenges_prompt_2 = """

#Challenges for the drug discovery program related to gamma secretase in Alzheimer's disease
                        
## Target Information 
### Develop a scientific rationale for the following:
                            
    **Given target:** gamma secretase
    **Given disease:** Alzheimer's disease
    **Given mode of action:** Inhibition of gamma secretase can reduce the production of amyloid-beta peptides and prevent the formation of amyloid plaques in the brain of Alzheimer's disease patients.
    
- Which patients would respond the therapy?
- Is the proposed mode of action on the target desirable and commercially viable in a clinical setting?
- What are advantages and disadvantages of different therapeutic modalities (antibodies, small molecules, antisense oligonucleotides, PROTACs, molecular glue, peptide macrocycles, and so on) for tackling the target?

- Alternative indications:
- Describe alternative indication for modulators of the target and explain why.


"""

task2_prompt = """ Task 2: Develop a target assessment strategy for gamma secretase in Alzeheimer's disease in maximal 500 words.
## Target Information 
                           
    **Given target:**  gamma secretase
    **Given disease:** Alzheimer's disease
    **Given mode of action:** Inhibition of gamma secretase can reduce amyloid beta production and aggregation


Outline a 1-year Target Assessment (TA) to Lead Identification (LI) plan. Describe High Level TA-LI plans.
- Make an emphasis on key inflection points that will inform the feasibility of the project. 
- Address status of in-vitro platforms, translational in vivo models (mechanistic models, not necessarily so called 'disease models')
  and describe what needs to be established. Elaborate on tractability and major challenges for advancement in a drug discovery portfolio.
- Discuss potential biomarkers and readouts for efficacy and target engagement.
"""
# this prompt is far too long - doubt it will lead to a good answer
safety_prompt = """Task 3: Safety assessment

## Target Information 
### Develop a scientific rationale for the following:
                           
    **Given target:**  gamma secretase
    **Given disease:** Alzheimer's disease
    **Given mode of action:** Inhibition of gamma secretase can reduce amyloid beta production and aggregation



- Does the target show bias towards expression in the desired organ (e.g. CNS)?
- Is it specifically expressed in the organ (e.g. brain)?
- Are there disease specific expression databases?
- Are there tissue-selective isoforms of the target?
- Are there condition-specific isoforms of the target?
- What regulates the alternative splicing that makes one isoform versus the other?
- How large is the expression of the target in the mouse model intended for in vivo tests?
- Is major phenotype reported in target knockouts and/or expression of rodent models?
- Are there published differences in expression between human and rodent models.
- What are the species differences that could be used to interpret rodent safety data on the target?
- What are the peripherial safety risks (oncogensis)?
- Can the modulation of the target promote tumor formation?
- Is there a way to assess on-target safety concerns?
- What are the safety concerns in case of exaggerated pharmacology?
- Will it disrupt cellular functions (e.g. endosomes, lysosomes, nuclear, mitochondrial) function with all its safety liability?
- How large is the risk for immunogenicity (related to biologics/antibody based approaches)?
- If the target is an enzyme, do polymorphisms in the human gene alter the protein enzyme activity?
 """

all_prompts = [working_hypothesis_prompt, clinical_target_prompt, Challenges_prompt_1, Challenges_prompt_2, task2_prompt, safety_prompt]



newline_char = '\n'
def run_with_timeout(func, timeout, *args, **kwargs):
    q = Queue()  # Create a Queue object for interprocess communication
    # Create a process to execute the given function, passing the Queue and other *args, **kwargs as parameters
    p = Process(target=func, args=(q, *args), kwargs=kwargs)
    p.start()
    # Wait for the process to complete or time out
    p.join(timeout)
    if p.is_alive():
        print(f"{datetime.now()} [INFO] Function {str(func)} execution timed out ({timeout}s), terminating process...")
        p.terminate()  # Terminate the process
        p.join()  # Ensure the process has been terminated
        result = None  # In case of a timeout, we do not have a result
    else:
        print(f"{datetime.now()} [INFO] Function {str(func)} completed successfully")
        result = q.get()  # Retrieve the result from the queue
    return result


def get_query_wrapper(q, question, answer):
    result = get_query(question, answer)
    q.put(result)

def get_query(question, answer):
    query_prompt = '''
I want to verify the content correctness of the given question, especially the last sentences.
Please summarize the content with the corresponding question.
This summarization will be used as a query to search with Pubmed central to find related knowledge or pages.
The query should be short but specific to promise retrieval of related knowledge or pages.
You can also use search syntax to make the query short and clear enough for the search engine to find relevant language data.
Try to make the query as relevant as possible to the last few sentences in the content.
**IMPORTANT**
Just output the query directly. DO NOT add additional explanations or introducement in the answer unless you are asked to.
'''
    query = openai.chat.completions.create(
        model="gpt-4o-mini",
        messages=[
            {
                "role": "system",
                "content": chat_prompt
            },
            {
                "role": "user",
                "content": f"##Question: {question}\n\n##Content: {answer}\n\n##Instruction: {query_prompt}"
            }
        ],
        temperature=0.5
    ).choices[0].message.content
    return query



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

def get_revise_answer_wrapper(q, question, answer, content):
    result = get_revise_answer(question, answer, content)
    q.put(result)

def get_revise_answer(question, answer, retrieved_data):
    revise_prompt = '''
I want to revise the answer according to retrieved related text of the question from Pubmed central. You will receive the summary, full text, relevance score to the query, and the full reference. If you use any of the given articles, you **must** reference using the name of the PDF throughout the text - e.g., [PDF_NAME].
**DO NOT REMOVE OR ALTER ANY CITATIONS ONLY ADD NEW ONES.** 
You need to check whether the answer is correct.
If you find some errors in the answer, revise the answer to make it better.
If you find some necessary details are ignored, add it to make the answer more plausible according to the related text.
If you find that a part of the answer is correct and does not require any additional details, maintain that part of the answer unchanged. Directly output the original content of that part without any modifications.
**IMPORTANT**
Try to keep the structure (multiple paragraphs with its subtitles) in the revised answer and make it more structural for understanding.
Split the paragraphs with `\n\n` characters.
For each piece of information you use from the retrieved data, cite the source at the end of the relevant sentence or paragraph.
Just output the revised answer directly. DO NOT add additional explanations or announcement in the revised answer unless you are asked to.
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



def generate_diff_html(text1, text2):
    diff = unified_diff(text1.splitlines(keepends=True),
                        text2.splitlines(keepends=True),
                        fromfile='text1', tofile='text2')
def RAG(question, draft_paragraphs):
    answer = ""
    all_retrieved_data = []
    for i, p in enumerate(draft_paragraphs):
        print("=" * 80)
        print(f"{datetime.now()} [INFO] Processing part {i + 1}/{len(draft_paragraphs)}...")
        answer += '\n\n' + p

        print(f"{datetime.now()} [INFO] Generating corresponding query...")
        res = run_with_timeout(get_query_wrapper, 90, question, answer)
        if not res:
            print(f"{datetime.now()} [INFO] Skipping subsequent steps...")
            continue
        else:
            query = res

        # Using a placeholder for newline character handling in the query
        newline_char = '\n'
        query_display = query.replace(newline_char, ' ')
        print(f">>> {i}/{len(draft_paragraphs)} Query: {query_display}")
        print(f"{datetime.now()} [INFO] Querying local database...")

        db_response = query_engine.query(query)

        if not db_response:
            print(f"{datetime.now()} [INFO] Skipping subsequent steps due to no response...")
            continue

        # Process the response and extract metadata
        retrieved_data = []
        for idx, node in enumerate(db_response.source_nodes):
            metadata = node.node.metadata
            source_id = f"SOURCE_{idx+1}"
            print(f"{datetime.now()} ############[INFO] Retrieved source {node.node.text} from source id {source_id} ########")
            retrieved_data.append({
                "source_id": source_id,
                "content": node.node.text,
                "metadata": metadata
            })
            
        all_retrieved_data.extend(retrieved_data)

        print(f"{datetime.now()} [INFO] Modifying answer based on database content...")
        res = run_with_timeout(get_revise_answer_wrapper, 90, question, answer, retrieved_data)
        if not res:
            print(f"{datetime.now()} [INFO] Skipping subsequent steps...")
        else:
            diff_html = generate_diff_html(answer, res)
            display(HTML(diff_html))
            answer = res
            print(f"{datetime.now()} [INFO] Answer modification completed")

    return answer, all_retrieved_data

def replace_source_ids_with_metadata(answer, retrieved_data):
    for source in retrieved_data:
        source_id = source['source_id']
        metadata = source['metadata']
        # Create a citation string from the metadata
        citation = f"(Source: {metadata.get('source', 'Unknown')}, Page: {metadata.get('page', 'N/A')})"
        # Replace the [SOURCE_ID] placeholder with the citation
        answer = answer.replace(f"[{source_id}]", citation)
    return answer


def split_draft(draft, split_char='\n\n'):
    # split_char: '\n\n'
    draft_paragraphs = draft.split(split_char)
    # print(f"The draft answer has {len(draft_paragraphs)}")
    return draft_paragraphs

def get_draft(question):
    # Getting the draft answer
    draft_prompt = '''
IMPORTANT:
Try to answer this question/instruction with step-by-step thoughts and make the answer more structural.
Use `\n\n` to split the answer into several paragraphs.
Just respond to the instruction directly. DO NOT add additional explanations or introducement in the answer unless you are asked to.
'''
    draft = openai.chat.completions.create(
        model="gpt-4o-mini",
        messages=[
            {
                "role": "system",
                "content": chat_prompt
            },
            {
                "role": "user",
                "content": f"{question}" + draft_prompt
            }
        ],
        temperature=0.5
    ).choices[0].message.content
    return draft

def get_draft_tot_inital(question, num_agents=3):

    draft_prompt = '''
IMPORTANT:
Try to answer this question/instruction with step-by-step thoughts and make the answer more structural.
Use `\n\n` to split the answer into several paragraphs.
Just respond to the instruction directly. DO NOT add additional explanations or introducement in the answer unless you are asked to.
'''

    refine_prompt = '''
Referencing the answers provided by all agents, synthesize a more detailed and comprehensive response by integrating all relevant details from these answers. Ensure logical coherence and provide ONLY THE MERGED ANSWER AS THE OUTPUT, omitting any discussion of the comparison process or analytical thoughts.
'''

    agents_drafts = []

    # Loop to generate different initial answers
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

        print(f"{datetime.now()} [INFO] Processing draft...")
        draft_paragraphs = split_draft(draft)
        print(f"{datetime.now()} [INFO] Draft split into {len(draft_paragraphs)} parts")

        # Modify using RAG
        draft_modified = RAG(question, draft_paragraphs)

        # Add each generated draft to the list
        agents_drafts.append(f"Agent{i+1}: {draft_modified}")

        print(f"{datetime.now()} [INFO] Agent{i + 1}/{num_agents} retrieved draft...")




    # Integrate and process previous answers
    agents_input = '\n\n'.join(agents_drafts) + '\n\n' + refine_prompt

    # Generate the integrated answer
    final_draft = openai.chat.completions.create(
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

    print(f"{datetime.now()} [INFO] Retrieved integrated draft...")

    return final_draft

def get_draft_tot(question, previous_answer, num_agents=3):
    # Update the draft answer prompt to include the question and previous answer
    draft_prompt = f'''
Base your response on the provided question and the previous answer. Expand the answer by adding more details to enhance its comprehensiveness. Ensure that the expansion maintains logical coherence and enriches the details, making the response more thorough and well-structured.
Question: {question}
Previous Answer: {previous_answer}
IMPORTANT:
Answer the full question with step-by-step thoughts and make the answer more structural.
Use `\n\n` to split the answer into several paragraphs.
Just respond to the instruction directly. DO NOT add additional explanations or introducement in the answer unless you are asked to.
'''

    refine_prompt = '''
Referencing the answers provided by all agents, synthesize a more detailed and comprehensive response by integrating all relevant details from these answers. Ensure logical coherence and provide ONLY THE MERGED ANSWER AS THE OUTPUT, omitting any discussion of the comparison process or analytical thoughts.
'''


    agents_drafts = []

    # Loop to generate initial different responses
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

        print(f"{datetime.now()} [INFO] Processing draft...")
        draft_paragraphs = split_draft(draft)
        print(f"{datetime.now()} [INFO] Draft split into {len(draft_paragraphs)} parts")

        # Modify using RAG
        draft_modified = RAG(question, draft_paragraphs)

        # Add each generated draft to the list
        agents_drafts.append(f"Agent{i + 1}: {draft_modified}")

        print(f"{datetime.now()} [INFO] Agent{i + 1}/{num_agents} retrieved draft...")

    # Integrate and process previous responses
    agents_input = '\n\n'.join(agents_drafts) + '\n\n' + refine_prompt

    # Generate the integrated answer
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

    print(f"{datetime.now()} [INFO] Retrieved integrated draft...")

    # Merge the integrated answer with the previous answer, prioritizing the previous answer with supplementary details from the new answer
    revise_prompt = f'''
Based on the original answer and an additional supplementary answer, generate a response that is richer in detail and logically coherent. Review the original answer:
1. If any part of the answer is correct and requires no further details, retain that portion unchanged and output it directly as it is.
2. For parts that may be improved or lack necessary details, enhance them by integrating information from the supplementary answer to make the response more comprehensive and accurate.
3. If you identify any errors within the answers, correct these errors while ensuring that the revised content remains logically coherent.
Original Answer: {previous_answer}
Supplementary Answer: {final_draft_raw}

**IMPORTANT**
Ensure the revised answer maintains a structured format (multiple paragraphs with subtitles) for better clarity. Separate the paragraphs with `\n\n` characters. Output only the enhanced answer directly, without any extra explanations or announcements unless specifically requested.
'''

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

    # Return the final merged draft
    return final_draft

@weave.op()
def ratt(question):
    step_num = num_steps
    print(f"{datetime.now()} [INFO] Retrieving Step 1 draft...")
    draft = get_draft_tot_inital(question, num_agents)
    print(f"{datetime.now()} [INFO] Step 1 draft returned")
    print(f"##################### DRAFT #######################")
    print(draft)
    print(f"#####################  END  #######################")

    print(f"{datetime.now()} [INFO] Processing draft...")
    draft_paragraphs = split_draft(draft)
    print(f"{datetime.now()} [INFO] Draft split into {len(draft_paragraphs)} parts")

    answer_first_state, retrieved_data = RAG(question, draft_paragraphs)
    answer_first_state = replace_source_ids_with_metadata(answer_first_state, retrieved_data)

    previous_answer = answer_first_state

    each_step_drafts = [f"Step 1 \n: {previous_answer}"]

    for iteration in range(1, step_num):
        print(f"{datetime.now()} [INFO] Retrieving Step {iteration + 1} draft...")
        draft = get_draft_tot(question, previous_answer, num_agents=num_agents)
        print(f"{datetime.now()} [INFO] Step {iteration + 1} draft returned")
        print(f"##################### DRAFT #######################")
        print(draft)
        print(f"#####################  END  #######################")

        print(f"{datetime.now()} [INFO] Processing draft...")
        draft_paragraphs = split_draft(draft)
        print(f"{datetime.now()} [INFO] Draft split into {len(draft_paragraphs)} parts")

        final_answer, retrieved_data = RAG(question, draft_paragraphs)
        final_answer = replace_source_ids_with_metadata(final_answer, retrieved_data)

        each_step_drafts.append(f"Step {iteration + 1} \n: {final_answer}")

        # Update previous_answer for the current iteration's response
        previous_answer = final_answer

    # Obtain the COT answer for baseline comparison
    draft_cot = get_draft(question)

    if final_output_mode == 'combine_each_step':
        final_draft = '\n\n'.join(each_step_drafts)
        refine_prompt = f'''
Referencing the answers provided by each step, synthesize a more detailed and comprehensive response by integrating all relevant details from these answers. Ensure logical coherence and provide ONLY THE MERGED ANSWER AS THE OUTPUT, omitting any discussion of the comparison process or analytical thoughts.
'''
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
            temperature=0.5
        ).choices[0].message.content

    return draft_cot, previous_answer

def replace_source_ids_with_metadata(answer, retrieved_data):
    for source in retrieved_data:
        source_id = source['source_id']
        metadata = source['metadata']
        # Create a citation string from the metadata
        citation = f"(Source: {metadata.get('source', 'Unknown')}, Page: {metadata.get('page', 'N/A')})"
        # Replace the [SOURCE_ID] placeholder with the citation
        answer = answer.replace(f"[{source_id}]", citation)
    return answer


def main():
    global query_engine
    if os.path.exists("./storage_enhanced"): # Check if the index has already been created
        print("Loading existing index...")
        storage_context = StorageContext.from_defaults(persist_dir="./storage_enhanced")
        sentence_index = load_index_from_storage(storage_context)
    else:
        print("Creating new index...")
        all_docs = []

            # Load all PDFs from the directory
        for filename in os.listdir(pdf_directory):
                if filename.endswith(".pdf"):
                    pdf_path = os.path.join(pdf_directory, filename)
                    
                    loader = PyPDFLoader(pdf_path)
                    all_pages = loader.load()
                                    
                    # Combine the text from the selected pages
                    doc_text = "\n\n".join([page.page_content for page in all_pages])
                    
                    # Create a Document object with metadata
                    doc = Document(text=doc_text, metadata={"source": filename})
                    print(f"Loaded PDF document: {filename}")
                    
                    all_docs.append(doc)

        print(f"Loaded {len(all_docs)} PDF documents.")

        # Set up the embedding model and Settings
        embedder = OpenAIEmbeddings(model="text-embedding-3-small")
        Settings.llm = OpenAI(model="gpt-4o-mini") 
        Settings.embed_model = embedder

            # Create the sentence window node parser with default settings
        node_parser = SentenceWindowNodeParser.from_defaults(
                window_size=3,
                window_metadata_key="window",
                original_text_metadata_key="original_text",
            )

            # Process all documents
        all_nodes = []
        for doc in all_docs:
                nodes = node_parser.get_nodes_from_documents([doc])
                all_nodes.extend(nodes)

            # Create the VectorStoreIndex
        sentence_index = VectorStoreIndex(all_nodes)

        print("VectorStoreIndex created successfully.")

        # Persist the index to disk
        sentence_index.storage_context.persist("./storage")

        print("Index persisted to disk.")

        # Now you can use sentence_index for querying, whether it was loaded or newly created
    query_engine = sentence_index.as_query_engine(
        similarity_top_k=4,
        # the target key defaults to `window` to match the node_parser's default
        node_postprocessors=[
            MetadataReplacementPostProcessor(target_metadata_key="window")
        ],
    )

    # # Printing out the values for demonstration
    print("Number of Agents:", num_agents)
    print("Number of Steps:", num_steps)
    print("Final Output Mode:", final_output_mode)

    weave.init(project_name="Ratt_Embedded_Enhanced")
    answer_cot, answer_ratt = ratt(task2_prompt)

    print(f"COT Answer:{answer_cot}")
    print(f"RATT Answer:{answer_ratt}")

    # window_response = query_engine.query(
    #     "Inhibition of gamma secretase for Alzheimer's disease treatment",
    # )
    # print(window_response)

    # for i, node in enumerate(window_response.source_nodes):
    #     print(f"Source {i + 1}:")
    #     print(node.node.metadata['window'])  # or node.node.text, depending on how it's stored
    #     print("\n")





if __name__ == '__main__':
    print("Script is running...")
    try:
        main()
    except Exception as e:
        print(f"An error occurred: {e}")
        import traceback
        traceback.print_exc()

