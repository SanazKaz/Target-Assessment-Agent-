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


gpt3_config = {
    "cache_seed": False,  # change the cache_seed for different trials
    "temperature": 0,
    "config_list": config_list,
    "timeout": 120,
}

initialiser = autogen.UserProxyAgent(
    name="Initialiser",
    human_input_mode="NEVER",
)

Moderator = autogen.AssistantAgent(
    name="Moderator",
    llm_config=config_list[0],
    human_input_mode="",
    system_message="""You are 'Moderator', overseeing the drug discovery process. 
                    ensure all aspects of the prompt are addressed, and synthesize a comprehensive final report. You will:
                    Your role is to coordinate the discussion between specialized agents,  ensuring all aspects of the prompt are addressed.

                    1. Analyze the given prompt and break it down into subtasks for each agent.
                    2. Assign tasks to appropriate agents and manage their interactions.
                    3. Ensure all required sections (scientific rationale, target assessment strategy, and safety assessment) are thoroughly covered.
                    4. Compile and summarize the inputs from all agents into a cohesive final report.
                    5. Ensure all claims are properly referenced and the report follows the specified format.
                    6. Present the final output, including a complete list of references."""
                    )

scientific_rational = autogen.AssistantAgent(
    name="SAR",
    llm_config=config_list[0],

    system_message="""You are 'SAR', an expert in target discovery and clinical target rationale. Your role is to develop a comprehensive scientific rationale for the given target in the specified disease. 

    
    You will:
                    Formulate a working hypothesis based on the provided information.
                    Analyze the clinical target rationale, focusing on human biology evidence.
                    Evaluate the relevance of the target to the disease and its physiological processes.
                    Assess genetic links, phenotypes, and clinical evidence related to the target.
                    Identify challenges in the drug discovery program for this target.
                    Suggest alternative indications for target modulators.
                    Use the use_pubmed_api function TWICE to find and cite relevant scientific literature.
                    Ensure your response covers all points requested in the prompt, with a minimum of 2000 words for the main rationale.""",)


safety_officer = autogen.AssistantAgent(
    name="SAFEA",
    system_message="""You are 'SAFEA' an expert in safety when it comes to target discovery. Your role is to conduct a thorough safety assessment of the proposed drug target. You will:

                    Evaluate target expression patterns and tissue specificity.
                    Analyze potential species differences between humans and animal models.
                    Assess on-target and off-target safety concerns, including peripheral risks.
                    Evaluate the risk of exaggerated pharmacology and immunogenicity.
                    Consider the impact of genetic polymorphisms on target function.
                    Provide a comprehensive safety profile addressing all points in the prompt's safety assessment section""",
    llm_config=config_list[0],
)
target_assessment = autogen.AssistantAgent(
    name="TAG",
    llm_config=config_list[0],
    system_message="""You are 'TAG' a specialist in assessing targets for drug discovery. Your role is to develop and outline a comprehensive 1-year Target Assessment (TA) to Lead Identification (LI) plan for the proposed drug target. Your responsibilities include:

                    Analyze and highlight key inflection points that will inform the project's feasibility.
                    Assess the status of in-vitro platforms and translational in vivo models.
                    Identify and describe what needs to be established in terms of platforms and models.
                    Elaborate on the target's tractability and major challenges for advancement in a drug discovery portfolio.
                    Discuss and propose potential biomarkers and readouts for efficacy and target engagement.
                    Outline a high-level TA-LI plan, focusing on critical milestones and decision points.
                    Evaluate the commercial viability and desirability of the proposed mode of action on the target in a clinical setting.
                    Assess the advantages and disadvantages of different therapeutic modalities for tackling the target.
                    Use the use_pubmed_api function to find relevant target assessment data and literature.
                    Provide a concise yet comprehensive target assessment strategy addressing all points in the prompt's target assessment section, within a 500-word limit.

                    Your assessment should be evidence-based, drawing from the latest research and industry best practices in target assessment. 
                    
                    Be prepared to interact with other specialists to ensure a well-rounded evaluation of the target""",
                    )

    # literature_agent = autogen.AssistantAgent(
    #     name="Literature_Agent",
    #     llm_config=config_list[0],
    #     system_message="""You provide relevant literature and references related to the prompt. 
    #                     You summarize the key points and provide a list of references.Use the use_pubmed_api function
    #                     to find relevant scientific literature when needed.""",)




agent_list=[scientific_rational, safety_officer, target_assessment]

    
groupchat = GroupChat(
    agents=[initialiser, Moderator, scientific_rational, safety_officer, target_assessment],
    messages=[],
    max_round=12,
    speaker_selection_method="round_robin",
    send_introductions=True,
)

manager = autogen.GroupChatManager(
    groupchat=groupchat,
    llm_config=gpt3_config,)


def use_pubmed_api(query: str , max_results=10):
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

# havent done it for moderator as it interferes with other specialist agents by suggesting tool call.
# A few options - 

# Register the use_pubmed_api function for all agents
for agent in agent_list:
    autogen.agentchat.register_function(
        f=use_pubmed_api,  
        caller=agent,
        executor=agent,
        name="use_pubmed_api",
        description="""
        Search PubMed for scientific literature related to the given query.

        This function should be used when you need to retrieve recent, peer-reviewed scientific 
        information from literature:
        - Finding evidence to support scientific claims
        - Gathering information on recent advancements in a specific area of biomedical research
        - Identifying key papers or authors in a particular field
        - Checking the current state of knowledge on a specific topic
        - You need up-to-date, scientifically validated information
        - You want to cite specific papers to support your arguments
        - You need to explore the current research landscape on a topic

        Returns: A formatted string containing details of the found papers, including titles, DOI, abstracts.
        """
    )




def process_message(message):
    # Clear the group chat history

    
    #groupchat.messages.clear()
    
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