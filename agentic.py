import os
from dotenv import load_dotenv
import json
import autogen
import tempfile
from autogen.coding import LocalCommandLineCodeExecutor
from autogen import GroupChat
from APIs.pubmed import PubMedAPI # Import the PubMedAPI class for literature search



# Load environment variables from .env file

load_dotenv()

openai_api_key = os.getenv("OPENAI_API_KEY")

# Load the config list from JSON
config_list = autogen.config_list_from_json(
    "OA_config_list.json",
    filter_dict={"model": ["gpt-3.5-turbo"]}
)

# Add the API key to each config in the list
for config in config_list:
    config["api_key"] = openai_api_key


pubmed_api =PubMedAPI(email="sanazkazemi@hotmail.com") # for enterez to contact you if neccessary

def use_pubmed_api(query, max_results=10):
    """
        Search PubMed for scientific literature related to the given query.

        This function should be used when you need to retrieve recent, peer-reviewed scientific 
        information from biomedical literature. It's particularly useful for:
        - Finding evidence to support scientific claims
        - Gathering information on recent advancements in a specific area of biomedical research
        - Identifying key papers or authors in a particular field
        - Checking the current state of knowledge on a specific topic

        Use this function when:
        - You need up-to-date, scientifically validated information
        - You want to cite specific papers to support your arguments
        - You need to explore the current research landscape on a topic

        Do not use this function when:
        - You need information from non-scientific sources
        - The query is not related to biomedical or life sciences
        - You require full-text articles (this API only provides abstracts)
        - You need information from a specific paper (use DOI or PMID instead)
        - You're looking for general knowledge that doesn't require scientific citation

        Args:
        query (str): The search query. Can include Boolean operators (AND, OR, NOT) and 
                    field tags (e.g., [Title], [Author], [Journal]).
        max_results (int, optional): Maximum number of results to return. Default is 5.
                                    Increasing this number will increase the API call duration.

        Returns:
        str: A formatted string containing details of the found papers, including titles, 
            authors, journal names, publication years, DOIs, and abstracts.

        Example:
        >>> result = use_pubmed_api("ACE inhibitors hypertension", 3)
        >>> print(result)
        Title: Comparative Effectiveness of ACE Inhibitors and ARBs in Hypertension Treatment
        Authors: Smith J, Johnson M, Williams R
        Journal: Journal of Hypertension, 2023
        DOI: 10.1000/jht.2023.1234
        Abstract: This study compares the effectiveness of ACE inhibitors and ARBs in...

        Note:
        - This function makes real-time API calls to PubMed. Use it judiciously to avoid 
        overwhelming the server or exceeding usage limits.
        - The results are based on PubMed's relevance ranking and may not always return 
        the most recent papers first.
        - Always critically evaluate the returned information and cross-reference when necessary.
        """
    
    print(f"PubMed API called with query: {query}, max_results: {max_results}")
    papers = pubmed_api.query(query, max_results)
    print(f"Retrieved {len(papers)} papers from PubMed")
    result = pubmed_api.format_results(papers)
    print("PubMed API call completed")
    
    return result
    

temp_dir = tempfile.TemporaryDirectory()
executor = LocalCommandLineCodeExecutor(
    timeout=10,  # Timeout for each code execution in seconds.
    work_dir=temp_dir.name,  # Use the temporary directory to store the code files.
)

gpt3_config = {
    "cache_seed": False,  # change the cache_seed for different trials
    "temperature": 0,
    "config_list": config_list,
    "timeout": 120,
}

initialiser = autogen.UserProxyAgent(
    name="Initialiser",
    human_input_mode="ALWAYS",
)

Moderator = autogen.AssistantAgent(
    name="Moderator",
    llm_config=config_list[0],
    human_input_mode="TERMINATE",
    system_message="""You are the Moderator overseeing the drug discovery process. Your role is to coordinate the discussion between specialized agents, ensure all aspects of the prompt are addressed, and synthesize a comprehensive final report. You will:

                    1. Analyze the given prompt and break it down into subtasks for each agent.
                    2. Assign tasks to appropriate agents and manage their interactions.
                    3. Ensure all required sections (scientific rationale, target assessment strategy, and safety assessment) are thoroughly covered.
                    4. Compile and summarize the inputs from all agents into a cohesive final report.
                    5. Ensure all claims are properly referenced and the report follows the specified format.
                    6. Present the final output, including a complete list of references."""
                    )

scientific_rational = autogen.AssistantAgent(
    name="Scientific_Rational",
    llm_config=config_list[0],
    system_message="""You are an expert in drug discovery and clinical target rationale. Your role is to develop a comprehensive scientific rationale for the given target in the specified disease. You will:
                    Formulate a working hypothesis based on the provided information.
                    Analyze the clinical target rationale, focusing on human biology evidence.
                    Evaluate the relevance of the target to the disease and its physiological processes.
                    Assess genetic links, phenotypes, and clinical evidence related to the target.
                    Identify challenges in the drug discovery program for this target.
                    Suggest alternative indications for target modulators.
                    Use the use_pubmed_api function to find and cite relevant scientific literature.
                    Ensure your response covers all points requested in the prompt, with a minimum of 2000 words for the main rationale.""",)


safety_officer = autogen.AssistantAgent(
    name="Safety_Assistant",
    system_message="""You are a safety expert in drug discovery. Your role is to conduct a thorough safety assessment of the proposed drug target. You will:

                    Evaluate target expression patterns and tissue specificity.
                    Analyze potential species differences between humans and animal models.
                    Assess on-target and off-target safety concerns, including peripheral risks.
                    Evaluate the risk of exaggerated pharmacology and immunogenicity.
                    Consider the impact of genetic polymorphisms on target function.
                    Use the use_pubmed_api function to find relevant safety data and literature.
                    Provide a comprehensive safety profile addressing all points in the prompt's safety assessment section.""",
    llm_config=config_list[0],
)
target_assessment = autogen.AssistantAgent(
    name="Target_Assessment",
    llm_config=config_list[0],
    system_message="""You are a safety expert in drug discovery. Your role is to conduct a thorough safety assessment of the proposed drug target. You will:
                    Evaluate target expression patterns and tissue specificity.
                    Analyze potential species differences between humans and animal models.
                    Assess on-target and off-target safety concerns, including peripheral risks.
                    Evaluate the risk of exaggerated pharmacology and immunogenicity.
                    Consider the impact of genetic polymorphisms on target function.
                    Use the use_pubmed_api function to find relevant safety data and literature.
                    Provide a comprehensive safety profile addressing all points in the prompt's safety assessment section.""",
                    )

    # literature_agent = autogen.AssistantAgent(
    #     name="Literature_Agent",
    #     llm_config=config_list[0],
    #     system_message="""You provide relevant literature and references related to the prompt. 
    #                     You summarize the key points and provide a list of references.Use the use_pubmed_api function
    #                     to find relevant scientific literature when needed.""",)

def state_transition(last_speaker, groupchat):
    print(f"checking if {last_speaker.name} used the pubmed api...")
    if hasattr(last_speaker, 'last_function_call') and last_speaker.last_function_call:
        if "use_pubmed_api" in last_speaker.last_function_call:
            print(f"{last_speaker.name} used the PubMed API")
            print(f"Query: {last_speaker.last_function_call['use_pubmed_api']['args'][0]}")
            print(f"Results: {last_speaker.last_function_call['use_pubmed_api']['return_value']}")
        else:
            print(f"{last_speaker.name} did not use the PubMed API")

    if last_speaker.name == "Initialiser":
        print("Next speaker: Moderator")
        return next(agent for agent in groupchat.agents if agent.name == "Moderator")
    elif last_speaker.name == "Moderator":
        print("Next speaker: Scientific_Rational")
        return next(agent for agent in groupchat.agents if agent.name == "Scientific_Rational")
    elif last_speaker.name == "Scientific_Rational":
        print("Next speaker: Safety_Assistant")
        return next(agent for agent in groupchat.agents if agent.name == "Safety_Assistant")
    elif last_speaker.name == "Safety_Assistant":
        print("Next speaker: Target_Assessment")
        return next(agent for agent in groupchat.agents if agent.name == "Target_Assessment")
    elif last_speaker.name == "Target_Assessment":
        print("Next speaker: Moderator")
        return next(agent for agent in groupchat.agents if agent.name == "Moderator")
    else:
        print(f"Unexpected last speaker: {last_speaker.name}. Ending conversation.")
        return None


agent_list=[initialiser, scientific_rational, safety_officer, target_assessment, Moderator]

for agent in agent_list:
    agent.register_function(
        function_map = {
            "use_pubmed_api": use_pubmed_api,
        })
    print(f"Registered use_pubmed_api function for agent: {agent.name}")
    
groupchat = GroupChat(
    agents=[initialiser,Moderator, scientific_rational, safety_officer, target_assessment],
    messages=[],
    max_round=6,
    speaker_selection_method=state_transition,
    send_introductions=True,
)


manager = autogen.GroupChatManager(
    groupchat=groupchat,
    llm_config=gpt3_config,
)

def process_message(message):
    # Clear the group chat history
    groupchat.messages.clear()
    
    # Initiate the chat using the initialiser
    initialiser.initiate_chat(
        manager,
        message=message,
        clear_history=True
    )
    
    # Extract the full conversation from the groupchat messages
    chat_history = groupchat.messages
    
    # Format the conversation
    formatted_conversation = []
    for msg in chat_history:
        if msg['content'].strip():  # Only add non-empty messages
            formatted_conversation.append({
                'speaker': msg['name'],
                'message': msg['content']
            })
    
    return formatted_conversation