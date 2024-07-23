summary_prompt = ("""
    Summarize the excerpt below to help answer a question.\n\n
    Excerpt from {citation}\n\n----\n\n{text}\n\n----\n\n
    Question: {question}\n\n
    Select papers that may help answer the question below. 
    Papers are listed as $KEY: $PAPER_INFO. 
    Return a list of keys, separated by commas. 
    'Return None, if no papers are applicable. '
    Choose papers that are relevant, from reputable sources, and timely 
    Do not directly answer the question, instead summarize to give evidence to help 
    Stay detailed; report specific numbers, equations, or 
    direct quotes (marked with quotation marks). Reply Not applicable if the
    excerpt is irrelevant. At the end of your response, provide an integer score 
    from 1-10 on a newline indicating relevance to question. Explain your score.
    'Respond with the following JSON format:

{{
  summary: ...,
  relevance_score: ...
}}
    \n\nRelevant Information Summary ({summary_length}):
)
                  
where `summary` is relevant information from text - {summary_length} words and `relevance_score` is the relevance of `summary` to answer question (out of 10).

""")
