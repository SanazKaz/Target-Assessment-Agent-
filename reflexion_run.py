import getpass
import os
from langchain.pydantic_v1 import BaseModel, Field
from langchain.tools import StructuredTool, tool
from langchain_openai import ChatOpenAI
from langchain_openai import OpenAIEmbeddings
from APIs.pubmed import Pubmed_API_langchain
import time
from langchain_core.messages import HumanMessage, ToolMessage
from langchain_core.output_parsers.openai_tools import PydanticToolsParser
from langchain_core.prompts import ChatPromptTemplate, MessagesPlaceholder
from langchain_core.pydantic_v1 import BaseModel, Field, ValidationError
import json
import datetime
from langchain_core.tools import StructuredTool
from langgraph.prebuilt import ToolNode
import random
from typing import Literal
from langgraph.graph.message import add_messages
from typing import Annotated
from langgraph.graph import END, MessageGraph, START
from typing_extensions import TypedDict
from langchain_text_splitters import SpacyTextSplitter
from typing import Literal
from langgraph.graph.message import add_messages
from typing import Annotated
from langgraph.graph import END, MessageGraph, START
from typing_extensions import TypedDict




working_hypothesis_prompt_GSEC = """ 
# Scientific Rationale for Gamma Secretase in Alzheimer's Disease


## Target Information 
### Develop a scientific rationale for the following:
                             
    **Given target:** Gamma secretase
    **Given disease:** Alzheimer's disease
    **Given mode of action:** Gamma secretase is a multi-subunit protease complex that cleaves type I transmembrane proteins, including the amyloid precursor protein (APP) leading to the generation of amyloid-beta (Aβ) peptides.

##Context:
Aβ is a family of secreted peptides generated from the sequential cleavages of the type 1 membrane protein APP by beta-secretase (BACE) and gamma-secretase (GSEC), respectively. 
BACE cleaves APP in the luminal domain, releasing the N-terminal soluble APPβ domain and leaving the C-terminal fragment, APP-CTF, which remains in the membrane. 
Subsequently, the APP-CTF is recruited to GSEC, a complex comprising four subunits, including PS, which harbors the active site. GSEC first cuts APP-CTF at the epsilon-cleavage site located close to the inner leaflet of the membrane. 
This cleavage event produces either Aβ48 or Aβ49 and the APP intracellular domain (AICD). The membrane-retained Aβ48 or Aβ49 is then further processed by GSEC in a continuous cascade of proteolytical events at every third of fourth amino acid, where the N-terminal product of each reaction becomes the substrate for the next GSEC cleavage event.
Accordingly, GSEC processes APP-CTF along two main product lines, Aβ49 → 46 → 43 → 40 → 37… and Aβ48 → 45 → 42 → 38…, respectively. During this processing cascade, Aβ43 and shorter Aβ peptides stochastically escape further processing by GSEC and are released into the extracellular space. 
As a result, Aβ peptides varying from 30 to 43 amino acids in length are secreted into the extracellular space. Among all secreted Aβ, Aβ40 is the most abundant in human CSF, followed by Aβ38, Aβ42, and Aβ37. In cognitively normal individuals, Aβ42 and Aβ43 represent a smaller portion of the total secreted Aβ.
These longer forms of Aβ seed the formation of Aβ-amyloid aggregates, a key step in the formation of amyloid plaques (Veugelen et al., 2016), as illustrated in Figure 1. Aβ42, which is produced in higher amounts than Aβ43, is the most abundant Aβ in amyloid plaques (Welander et al., 2009).

## Task 1: Develop Scientific Rationale

### Working Hypothesis
- Detailed description of the idea
- Unmet medical need
- Suitability for combination therapy
- Predictive biomarkers
- Clinical relevance of existing biomarkers


"""

clinical_target_prompt_GSEC = """ #Clinical Target Rationale for Gamma Secretase in Alzheimer's Disease


## Target Information 
### Develop a scientific rationale for the following:
                             
    **Given target:** Gamma secretase
    **Given disease:** Alzheimer's disease
    **Given mode of action:** Gamma secretase is a multi-subunit protease complex that cleaves type I transmembrane proteins, including the amyloid precursor protein (APP) leading to the generation of amyloid-beta (Aβ) peptides.

##Context:
Aβ is a family of secreted peptides generated from the sequential cleavages of the type 1 membrane protein APP by beta-secretase (BACE) and gamma-secretase (GSEC), respectively. 
BACE cleaves APP in the luminal domain, releasing the N-terminal soluble APPβ domain and leaving the C-terminal fragment, APP-CTF, which remains in the membrane. 
Subsequently, the APP-CTF is recruited to GSEC, a complex comprising four subunits, including PS, which harbors the active site. GSEC first cuts APP-CTF at the epsilon-cleavage site located close to the inner leaflet of the membrane. 
This cleavage event produces either Aβ48 or Aβ49 and the APP intracellular domain (AICD). The membrane-retained Aβ48 or Aβ49 is then further processed by GSEC in a continuous cascade of proteolytical events at every third of fourth amino acid, where the N-terminal product of each reaction becomes the substrate for the next GSEC cleavage event.
Accordingly, GSEC processes APP-CTF along two main product lines, Aβ49 → 46 → 43 → 40 → 37… and Aβ48 → 45 → 42 → 38…, respectively. During this processing cascade, Aβ43 and shorter Aβ peptides stochastically escape further processing by GSEC and are released into the extracellular space. 
As a result, Aβ peptides varying from 30 to 43 amino acids in length are secreted into the extracellular space. Among all secreted Aβ, Aβ40 is the most abundant in human CSF, followed by Aβ38, Aβ42, and Aβ37. In cognitively normal individuals, Aβ42 and Aβ43 represent a smaller portion of the total secreted Aβ.
These longer forms of Aβ seed the formation of Aβ-amyloid aggregates, a key step in the formation of amyloid plaques (Veugelen et al., 2016), as illustrated in Figure 1. Aβ42, which is produced in higher amounts than Aβ43, is the most abundant Aβ in amyloid plaques (Welander et al., 2009).
 

 ### Clinical target rationale:
    - How relevant is the target location to the disease biology?
    - How is the target expression altered in human disease?
    - How is the target involved in the physiological process relevant to the disease?
    - Which phenotypes and genotypes were identified for the target?
    - How is the genetic link between the target and the disease?
    - Describe the evidence provided in clinics or by tools acting on the pathway where the target is involved.
    - Which kind of target modulation is required to treat the disease? """

Challenges_prompt_GSEC = """ #Challenges for the drug discovery program related to Gamma Secretase as a target in Alzheimer's Disease
    **Given target:** Gamma secretase
    **Given disease:** Alzheimer's disease
    **Given mode of action:** Gamma secretase is a multi-subunit protease complex that cleaves type I transmembrane proteins, including the amyloid precursor protein (APP) leading to the generation of amyloid-beta (Aβ) peptides.

##Context:
Aβ is a family of secreted peptides generated from the sequential cleavages of the type 1 membrane protein APP by beta-secretase (BACE) and gamma-secretase (GSEC), respectively. 
BACE cleaves APP in the luminal domain, releasing the N-terminal soluble APPβ domain and leaving the C-terminal fragment, APP-CTF, which remains in the membrane. 
Subsequently, the APP-CTF is recruited to GSEC, a complex comprising four subunits, including PS, which harbors the active site. GSEC first cuts APP-CTF at the epsilon-cleavage site located close to the inner leaflet of the membrane. 
This cleavage event produces either Aβ48 or Aβ49 and the APP intracellular domain (AICD). The membrane-retained Aβ48 or Aβ49 is then further processed by GSEC in a continuous cascade of proteolytical events at every third of fourth amino acid, where the N-terminal product of each reaction becomes the substrate for the next GSEC cleavage event.
Accordingly, GSEC processes APP-CTF along two main product lines, Aβ49 → 46 → 43 → 40 → 37… and Aβ48 → 45 → 42 → 38…, respectively. During this processing cascade, Aβ43 and shorter Aβ peptides stochastically escape further processing by GSEC and are released into the extracellular space. 
As a result, Aβ peptides varying from 30 to 43 amino acids in length are secreted into the extracellular space. Among all secreted Aβ, Aβ40 is the most abundant in human CSF, followed by Aβ38, Aβ42, and Aβ37. In cognitively normal individuals, Aβ42 and Aβ43 represent a smaller portion of the total secreted Aβ.
These longer forms of Aβ seed the formation of Aβ-amyloid aggregates, a key step in the formation of amyloid plaques (Veugelen et al., 2016), as illustrated in Figure 1. Aβ42, which is produced in higher amounts than Aβ43, is the most abundant Aβ in amyloid plaques (Welander et al., 2009).

### Challenges:
- Check the following idea for details on small molecule compounds: Developing small molecule modulators or inhibitors of gamma secretase for Alzheimer's disease treatment.
- Is a 'information driven approach' (IDA) strategy based on available small molecules possible?
- Which small molecular modulators of the target known?
- Which inhibitors, antagonists, agonists, negative allosteric modulators (NAM), positive allosteric modulators (PAM) are required for target modulation in the given disease? 
- Which patients would respond the therapy?
- Is the proposed mode of action on the target desirable and commercially viable in a clinical setting?
- What are advantages and disadvantages of different therapeutic modalities (antibodies, small molecules, antisense oligonucleotides, PROTACs, molecular glue, peptide macrocycles, and so on) for tackling the target?

- Alternative indications:
- Describe alternative indication for modulators of the target and explain why.


"""



working_hypothesis_prompt_BACE = """ 
# Scientific Rationale for BACE in Alzheimer's Disease


## Target Information 
### Develop a scientific rationale for the following:
                             
    **Given target:** β-Site APP cleaving enzyme 1 (BACE1)
    **Given disease:** Alzheimer's disease
    **Given mode of action:** BACE1 inhibition leads to reduced production of amyloid-β (Aβ) peptides, which are the main component of amyloid plaques in the brains of Alzheimer's patients.

## Task 1: Develop Scientific Rationale

### Working Hypothesis
- Detailed description of the idea
- Unmet medical need
- Suitability for combination therapy
- Predictive biomarkers
- Clinical relevance of existing biomarkers


"""

clinical_target_prompt_BACE = """ #Clinical Target Rationale for BACE in Alzheimer's Disease


## Target Information 
### Develop a scientific rationale for the following:
                           
    **Given target:** β-Site APP cleaving enzyme 1 (BACE1)
    **Given disease:** Alzheimer's disease
    **Given mode of action:** BACE1 inhibition leads to reduced production of amyloid-β (Aβ) peptides, which are the main component of amyloid plaques in the brains of Alzheimer's patients.

##Context:
Aβ is a family of secreted peptides generated from the sequential cleavages of the type 1 membrane protein APP by beta-secretase (BACE) and gamma-secretase (GSEC), respectively. 
BACE cleaves APP in the luminal domain, releasing the N-terminal soluble APPβ domain and leaving the C-terminal fragment, APP-CTF, which remains in the membrane. 
Subsequently, the APP-CTF is recruited to GSEC, a complex comprising four subunits, including PS, which harbors the active site. GSEC first cuts APP-CTF at the epsilon-cleavage site located close to the inner leaflet of the membrane. 
This cleavage event produces either Aβ48 or Aβ49 and the APP intracellular domain (AICD). The membrane-retained Aβ48 or Aβ49 is then further processed by GSEC in a continuous cascade of proteolytical events at every third of fourth amino acid, where the N-terminal product of each reaction becomes the substrate for the next GSEC cleavage event.
Accordingly, GSEC processes APP-CTF along two main product lines, Aβ49 → 46 → 43 → 40 → 37… and Aβ48 → 45 → 42 → 38…, respectively. During this processing cascade, Aβ43 and shorter Aβ peptides stochastically escape further processing by GSEC and are released into the extracellular space. 
As a result, Aβ peptides varying from 30 to 43 amino acids in length are secreted into the extracellular space. Among all secreted Aβ, Aβ40 is the most abundant in human CSF, followed by Aβ38, Aβ42, and Aβ37. In cognitively normal individuals, Aβ42 and Aβ43 represent a smaller portion of the total secreted Aβ.
These longer forms of Aβ seed the formation of Aβ-amyloid aggregates, a key step in the formation of amyloid plaques (Veugelen et al., 2016), as illustrated in Figure 1. Aβ42, which is produced in higher amounts than Aβ43, is the most abundant Aβ in amyloid plaques (Welander et al., 2009).
 

 ### Clinical target rationale:
    - How relevant is the target location to the disease biology?
    - How is the target expression altered in human disease?
    - How is the target involved in the physiological process relevant to the disease?
    - Which phenotypes and genotypes were identified for the target?
    - How is the genetic link between the target and the disease?
    - Describe the evidence provided in clinics or by tools acting on the pathway where the target is involved.
    - Which kind of target modulation is required to treat the disease? """

Challenges_prompt_BACE = """ #Challenges for the drug discovery program related to BACE as a target in Alzheimer's Disease
                           
    **Given target:** β-Site APP cleaving enzyme 1 (BACE1)
    **Given disease:** Alzheimer's disease
    **Given mode of action:** BACE1 inhibition leads to reduced production of amyloid-β (Aβ) peptides, which are the main component of amyloid plaques in the brains of Alzheimer's patients.

##Context:
Aβ is a family of secreted peptides generated from the sequential cleavages of the type 1 membrane protein APP by beta-secretase (BACE) and gamma-secretase (GSEC), respectively. 
BACE cleaves APP in the luminal domain, releasing the N-terminal soluble APPβ domain and leaving the C-terminal fragment, APP-CTF, which remains in the membrane. 
Subsequently, the APP-CTF is recruited to GSEC, a complex comprising four subunits, including PS, which harbors the active site. GSEC first cuts APP-CTF at the epsilon-cleavage site located close to the inner leaflet of the membrane. 
This cleavage event produces either Aβ48 or Aβ49 and the APP intracellular domain (AICD). The membrane-retained Aβ48 or Aβ49 is then further processed by GSEC in a continuous cascade of proteolytical events at every third of fourth amino acid, where the N-terminal product of each reaction becomes the substrate for the next GSEC cleavage event.
Accordingly, GSEC processes APP-CTF along two main product lines, Aβ49 → 46 → 43 → 40 → 37… and Aβ48 → 45 → 42 → 38…, respectively. During this processing cascade, Aβ43 and shorter Aβ peptides stochastically escape further processing by GSEC and are released into the extracellular space. 
As a result, Aβ peptides varying from 30 to 43 amino acids in length are secreted into the extracellular space. Among all secreted Aβ, Aβ40 is the most abundant in human CSF, followed by Aβ38, Aβ42, and Aβ37. In cognitively normal individuals, Aβ42 and Aβ43 represent a smaller portion of the total secreted Aβ.
These longer forms of Aβ seed the formation of Aβ-amyloid aggregates, a key step in the formation of amyloid plaques (Veugelen et al., 2016), as illustrated in Figure 1. Aβ42, which is produced in higher amounts than Aβ43, is the most abundant Aβ in amyloid plaques (Welander et al., 2009).

### Challenges:
- Check the following idea for details on small molecule compounds: Developing small molecule modulators or inhibitors of BACE for Alzheimer's disease treatment.
- Is a 'information driven approach' (IDA) strategy based on available small molecules possible?
- Which small molecular modulators of the target known?
- Which inhibitors, antagonists, agonists, negative allosteric modulators (NAM), positive allosteric modulators (PAM) are required for target modulation in the given disease? 
- Which patients would respond the therapy?
- Is the proposed mode of action on the target desirable and commercially viable in a clinical setting?
- What are advantages and disadvantages of different therapeutic modalities (antibodies, small molecules, antisense oligonucleotides, PROTACs, molecular glue, peptide macrocycles, and so on) for tackling the target?

- Alternative indications:
- Describe alternative indication for modulators of the target and explain why.


"""



working_hypothesis_prompt_ALPHA = """ 
# Scientific Rationale for 5-alpha reductase  in Huntington's Disease


## Target Information 
### Develop a scientific rationale for the following:
                             
    **Given target:**  5-alpha reductase
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

clinical_target_prompt_ALPHA = """ #Clinical Target Rationale for 5-alpha reductase in Huntington's Disease


## Target Information 
### Develop a scientific rationale for the following:
                           
    **Given target:** 5-alpha reductase
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

Challenges_prompt_ALPHA = """ 
#Challenges for the drug discovery program related to 5-alpha reductase  in Huntington's disease
                           
## Target Information 
### Develop a scientific rationale for the following:
                           
    **Given target:**  5-alpha reductase
    **Given disease:** Huntington's disease
    **Given mode of action:** Inhibition of 5-alpha reductase can promote neuronal survival and reduce mutant huntingtin protein aggregation


### Challenges:
- Check the following idea for details on small molecule compounds: Developing small molecule modulators or inhibitors of 5-alpha reductase  for Huntington's Disease treatment.
- Is a 'information driven approach' (IDA) strategy based on available small molecules possible?
- Which small molecular modulators of the target known?
- Which inhibitors, antagonists, agonists, negative allosteric modulators (NAM), positive allosteric modulators (PAM) are required for target modulation in the given disease? 

- Which patients would respond the therapy?
- Is the proposed mode of action on the target desirable and commercially viable in a clinical setting?
- What are advantages and disadvantages of different therapeutic modalities (antibodies, small molecules, antisense oligonucleotides, PROTACs, molecular glue, peptide macrocycles, and so on) for tackling the target?

- Alternative indications:
- Describe alternative indication for modulators of the target and explain why.


"""




llm = ChatOpenAI(model="gpt-4o", temperature=0)
embedding_model = OpenAIEmbeddings()

def _set_if_undefined(var: str) -> None:
    if os.environ.get(var):
        return
    os.environ[var] = getpass.getpass(var)


# Optional: Configure tracing to visualize and debug the agent
_set_if_undefined("LANGCHAIN_API_KEY")
os.environ["LANGCHAIN_TRACING_V2"] = "true"
os.environ["LANGCHAIN_PROJECT"] = "Reflexion_1"

_set_if_undefined("OPENAI_API_KEY")



@tool
def get_pubmed_papers_for_llm(query, max_results=4, email="your_email@example.com", sleep_time=3, used_papers=set()):
    """
    Retrieve papers from PubMed Central and format them for LLM processing, skipping previously used papers.
    
    Args:
    query (str): The search query for PubMed Central.
    max_results (int): Maximum number of papers to retrieve. Default is 4.
    email (str): Your email address for NCBI's records. Replace with your actual email.
    sleep_time (float): Time to sleep between API requests in seconds. Default is 3.
    used_papers (set): Set of paper IDs that have already been used and should be skipped.
    
    Returns:
    list: A list of dictionaries, each containing a paper's content and metadata.
    """
    print("Initializing PMC API")
    pmc_api = Pubmed_API_langchain(email=email)
    
    print(f"Searching PMC for: '{query}'")
    papers = pmc_api.query(query, max_results)
    print(f"Retrieved {len(papers)} papers")
    
    processed_papers = []
    
    for i, paper in enumerate(papers):
        if len(processed_papers) >= max_results:
            break
        
        if paper['pmc_id'] in used_papers:
            print(f"Skipping already used paper: {paper['pmc_id']}")
            continue
        
        if i > 0:
            time.sleep(sleep_time)
        
        full_content = f"Title: {paper['title']}\n\nAbstract: {paper['abstract']}\n\nFull Text: {paper['full_text']}"

        text_splitter = SpacyTextSplitter()
        tokens = text_splitter.split_text(full_content)
        token_count = len(tokens)
        
        if token_count > 100000:
            print(f"Truncating paper {paper['pmc_id']} from {token_count} to 100,000 tokens")
            tokens = tokens[:100000]
            full_content = ' '.join(tokens)
            token_count = 100000
            print(f"Truncated content: {full_content}")
        
        
        processed_paper = {
            'content': full_content,
            'metadata': {
                'paper_id': paper['pmc_id'],
                'authors': paper['authors'],
                'year': paper['year'],
                'journal': paper['journal']
            }
        }

        processed_papers.append(processed_paper)
        used_papers.add(paper['pmc_id'])

    # Save used papers
    with open("used_papers_reflexion.txt", "w") as f:
        f.write("\n".join(used_papers))
    
    # Print processing summary
    if processed_papers:
        print("\nPaper structure:")
        print(f"Content length: {len(processed_papers[-1]['content'])}")
        print(f"Metadata keys: {', '.join(processed_papers[-1]['metadata'].keys())}")

    else:
        print("No new papers were processed.")

    return processed_papers



class Reflection(BaseModel):
    missing: str = Field(description="Critique of what is missing.")
    superfluous: str = Field(description="Critique of what is superfluous and could be removed or improved with depth.")


class AnswerQuestion(BaseModel):
    """Provide an answer, each paragraph must have a reference from search results, reflection, and then follow up with search queries to improve the answer. Provide a detailed summary of relevant search"""
    answer: str = Field(description=""""You are a AI assistant with a background in drug discovery and specialise in writing scientific rationales for given targets and diseases, you cite your work based on the retrieved pubmed results at the end of relevant sentence as follows:
                        [Harold et al, Pubmed ID].""")
    
    reflection: Reflection = Field(description="Your reflection on the initial answer.")

    # best_query: str = Field(description="The best selected search query")
    search_queries: list[str] = Field(
        description="""Generate ONE diverse search QUERY using advanced search operators to
          improve the answer and bring back relevant results - this will be used in pubmed central so make it as relevant as possible 
        ##use this as an example - "BACE1" AND "Alzheimer's disease" AND ("Notch signaling" OR "GSEC development" OR "drug discovery").
        ###The search query should be related to different aspects of the prompt """),
    
    search_results: str = Field(description="A summary of relevant search_results with quantitative details and examples for use in the answer, add in metadata.")
class AnswerQuestion(BaseModel):
    answer: str = Field(description="Your answer here")
    reflection: Reflection = Field(description="Your reflection on the initial answer")
    search_queries: list[str] = Field(description="Generate ONE diverse search QUERY")
    search_results: str = Field(description="A summary of relevant search results with details and examples for use in the answer")

    def to_json(self):
        return json.dumps(self.dict(), ensure_ascii=False, indent=2)
class ResponderWithRetries:
    def __init__(self, runnable, validator):
        self.runnable = runnable
        self.validator = validator

    def respond(self, state: list):
        response = []
        for attempt in range(3):
            response = self.runnable.invoke(
                {"messages": state}, {"tags": [f"attempt:{attempt}"]}
            )
            try:
                self.validator.invoke(response)
                return response
            except ValidationError as e:
                state = state + [
                    response,
                    ToolMessage(
                        content=f"{repr(e)}\n\nPay close attention to the function schema.\n\n"
                        + self.validator.schema_json()
                        + " Respond by fixing all validation errors.",
                        tool_call_id=response.tool_calls[0]["id"],
                    ),
                ]
        return response

    def format_error_message(self, e):
        # Implement this method to format the error message
        return f"Validation error: {str(e)}"

    def format_error_message(self, error: ValidationError):
        # Extract relevant information from the error
        error_details = error.errors()
        formatted_errors = []
        for detail in error_details:
            field = detail['loc'][0] if detail['loc'] else 'unknown field'
            message = detail['msg']
            formatted_errors.append(f"- {field}: {message}")
        
        error_message = "Your response didn't meet the expected format. Please fix the following issues:\n"
        error_message += "\n".join(formatted_errors)
        error_message += "\n\nPlease provide your response in valid MARKDOWN format with the following structure:\n"
        error_message += '{\n  "answer": "Your Markdown-formatted answer here",\n  "references": ["Reference 1", "Reference 2", ...]\n}'
        
        return error_message
    
    import datetime

actor_prompt_template = ChatPromptTemplate.from_messages(
    [
        (
            "system",
            """You are an AI assistant specializing in disease and drug discovery research. Your task is to analyze scientific papers and provide a comprehensive response.
1. You will receive a scientific paper in this format:
- content: The paper  content. Use this to improve your answer each time.
- metadata: Paper ID, authors, year, and journal.
When using information from metadata, cite them as [Paper ID, Author] at the end of the relevant sentence. You cite every relevant sentence otherwise it will be considered plagiarism.
2. Suggest search queries for further research, using advanced techniques.
Example query: "BACE1" AND "Alzheimer's disease" AND ("Notch signaling" OR "GSEC development" OR "drug discovery")
Maintain a critical perspective and provide accurate, up-to-date information. Continuously refine your answer with each iteration.
3. Your response should be detailed, structured, and supported by references with a full reference list at the end.""",

        ),
        MessagesPlaceholder(variable_name="messages"),
        (
            "user",
            "\n\n<system>Reflect on the user's original question and the"
            " actions taken thus far. Respond using the {function_name} function.</reminder>",
        ),
    ]
).partial(
    time=lambda: datetime.datetime.now().isoformat(),
)
initial_answer_chain = actor_prompt_template.partial(
    first_instruction="Provide an extremely detailed and comprehensive structured answer to the question.",
    function_name=AnswerQuestion.__name__,
) | llm.bind_tools(tools=[AnswerQuestion])
validator = PydanticToolsParser(tools=[AnswerQuestion])

first_responder = ResponderWithRetries(
    runnable=initial_answer_chain, validator=validator
)

initial = first_responder.respond([HumanMessage(content=Challenges_prompt_ALPHA)])


from typing import Optional


revise_instructions = """ you are an AI assistant specializing in scientific literature analysis and synthesis, with a focus on biomedical research, particularly in Alzheimer's disease and drug discovery. Your task is to revise and improve your previous answer based on new information. Follow these steps:

1. Carefully review the new scientific papers provided in the search results.
2. Compare this new information with your previous answer.
3. Integrate the new information to enhance and expand your answer. Be sure to:
   - Update any outdated information
   - Add new relevant details, especially about current research and clinical trials
   - Strengthen your arguments with additional evidence
   - Address any gaps or weaknesses identified in your previous reflection

4. Maintain a clear and logical structure in your revised answer.
5. Highlight new insights or connections you've made based on the new information.
6. If you encounter conflicting information, discuss these discrepancies and provide a balanced view.
7. Ensure all major claims are supported by at least one reference, citing sources using Harvard Style.
8. After revising, reflect on how the new information has changed or enhanced your understanding.
9. Identify any remaining gaps in your knowledge or areas of uncertainty.
10. Suggest ONE search query using search notation that could address these gaps or uncertainties in the next iteration - be flexible by using OR but specific using shrot keywords with quotes.

Your goal is to continuously refine and improve your answer, making it more comprehensive, accurate, and up-to-date with each revision. 
Remember to maintain a critical and analytical perspective throughout this process. """

# Extend the initial answer schema to include references.
# Forcing citation in the model encourages grounded responses
class ReviseAnswer(AnswerQuestion):
    search_results: str = Field(description="""A detailed summary of **RELEVANT** information retrieved from the pubmed search and its metadata -
                                 give examples of details that are important to the original prompt.""")

    references: list[str] = Field(
        default_factory=list, # will always have an empty list if not provided by model
        description="Citations from search_results motivating your updated answer."
    )
    """Revise your original answer to the question. Provide an answer, reflection,

    cite your reflection with references in the Harvard style, and finally
    suggest 10 search queries based on your reflection to improve the answer - **the first query should be DIFFERENT every time**"""


initial_answer_chain = actor_prompt_template.partial(
    first_instruction="Provide an extremely detailed and comprehensive structured answer to the question.",
    function_name=AnswerQuestion.__name__,
) | llm.bind_tools(tools=[AnswerQuestion])
validator = PydanticToolsParser(tools=[AnswerQuestion])


revision_chain = actor_prompt_template.partial(
    first_instruction=revise_instructions,
    function_name=ReviseAnswer.__name__,
) | llm.bind_tools(tools=[ReviseAnswer])
revision_validator = PydanticToolsParser(tools=[ReviseAnswer])

revisor = ResponderWithRetries(runnable=revision_chain, validator=revision_validator)

import json

search_queries = []

revised = revisor.respond(
    [
        HumanMessage(content=Challenges_prompt_ALPHA),
        initial,
        ToolMessage(
            tool_call_id=initial.tool_calls[0]["id"],
            content=json.dumps(
                get_pubmed_papers_for_llm.invoke(
                    {"query": initial.tool_calls[0]["args"]["search_queries"][0]}
                )
            ),
        ),
    ]
)
revised

print(json.dumps(initial.tool_calls[0]["args"], indent=2))


from langchain_core.tools import StructuredTool
from langgraph.prebuilt import ToolNode
import random



used_queries = set() # Keep track of used queries to avoid duplicates

def run_queries(search_queries: list[str], **kwargs):
    """Run an unused query if available, otherwise use a random one."""
    global used_queries
    unused_queries = [q for q in search_queries if q not in used_queries]
    if unused_queries:
        selected_query = random.choice(unused_queries)
    elif search_queries:
        selected_query = random.choice(search_queries)
    else:
        return "No search queries provided."
    
    used_queries.add(selected_query)
    return get_pubmed_papers_for_llm.invoke({"query": selected_query})



tool_node = ToolNode(
    [
        StructuredTool.from_function(run_queries, name=AnswerQuestion.__name__),
        StructuredTool.from_function(run_queries, name=ReviseAnswer.__name__),
    ]
)


MAX_ITERATIONS = 7
builder = MessageGraph()

class State(TypedDict):
    messages: Annotated[list, add_messages]

builder.add_node("draft", first_responder.respond)


builder.add_node("execute_tools", tool_node)
builder.add_node("revise", revisor.respond)
# draft -> execute_tools
builder.add_edge("draft", "execute_tools")
# execute_tools -> revise
builder.add_edge("execute_tools", "revise")

# Define looping logic:


def _get_num_iterations(state: list):
    i = 0
    for m in state[::-1]:
        if m.type not in {"tool", "ai"}:
            break
        i += 1
    return i


def event_loop(state: list) -> Literal["execute_tools", "__end__"]:
    # in our case, we'll just stop after N plans
    num_iterations = _get_num_iterations(state)
    if num_iterations > MAX_ITERATIONS:
        return END
    return "execute_tools"


# revise -> execute_tools OR end
builder.add_conditional_edges("revise", event_loop)
builder.add_edge(START, "draft")
graph = builder.compile()



events = graph.stream(
    [HumanMessage(content=Challenges_prompt_ALPHA)],
    stream_mode="values",
)
for i, step in enumerate(events):
    print(f"Step {i}")
    step[-1].pretty_print()
