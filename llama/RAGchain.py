from langchain_community.embeddings import OllamaEmbeddings
from langchain_community.llms import Ollama
from langchain_community.vectorstores import FAISS
from langchain_community.document_loaders import PyPDFLoader
from langchain_text_splitters import RecursiveCharacterTextSplitter
from langchain_core.prompts import ChatPromptTemplate
from langchain_core.output_parsers import StrOutputParser
from langchain.chains import create_history_aware_retriever
from langchain_core.prompts import MessagesPlaceholder
from langchain_core.messages import HumanMessage, AIMessage



# Load your PDF documents
loader = PyPDFLoader("/Users/sanazkazeminia/Documents/LLM_Agent/science-in-the-age-of-ai-report.pdf")
document = loader.load()

# Create embeddings for the document chunks
embeddings = OllamaEmbeddings(model="llama3")

# Split the documents into chunks
text_splitter = RecursiveCharacterTextSplitter()
documents = text_splitter.split_documents(document) # Split the split documents into chunks


# Create a vector store from the document chunks and embeddings
vector = FAISS.from_documents(documents, embeddings)

# Initialize the Llama3 model
llm = Ollama(model="llama3")

output_parser = StrOutputParser()

# Create a chain of LLMs
prompt = ChatPromptTemplate.from_messages([
    MessagesPlaceholder(variable_name="chat_history"),
    ("system", "You are an expert in the field of AI for Science. I have a question for you."),
    ("user", "{input}"),
    ("user", "Given the above conversation, generate a search query to look up to get information relevant to the conversation")
])
retriever_chain = create_history_aware_retriever(llm, retriever, prompt)

chat_history = [HumanMessage(content="Can LangSmith help test my LLM applications?"), AIMessage(content="Yes!")]
retriever_chain.invoke({
    "chat_history": chat_history,
    "input": "Tell me how"
})

document_chain = create_stuff_documents_chain(llm, prompt)

retrieval_chain = create_retrieval_chain(retriever_chain, document_chain)


chat_history = [HumanMessage(content="Can LangSmith help test my LLM applications?"), AIMessage(content="Yes!")]
retrieval_chain.invoke({
    "chat_history": chat_history,
    "input": "Tell me how"
})