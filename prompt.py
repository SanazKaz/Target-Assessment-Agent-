def prompt_to_use(target, disease, mode_of_action, context, idea):

    prompt = f"""You are a AI assistant with a background in drug discovery.

    Given target: {target}
    Given disease: {disease}
    Given mode of action: {mode_of_action}

    Context:
    {context}

    Task 1: Develop a scientific rationale for {target} in {disease}.

    Highlight the working hypothesis for the clinical target rationale and human biology evidence by minimum 2000 words.

    Describe as much as possible the evidence in humans or in human tissue that link the target, target space or approach to the pathogenesis of interest.
    If known, also describe here the wanted mode of action with regards to desired clinical outcome.
    Please avoid including only pre-clinical data in this section.

    Use the following structure and provide a detailed description for each point:
    - Working hypothesis:
    - Create a detailed description of the following idea: {idea}
    - Is there are significant unmet medical need?
    - Is it suitable for combination therapy?
    - Which predictive biomarkers exist for the target related to the disease?
        - Provide a detailed description of existing clinical relevant biomarkers.

    - Clinical target rationale:
    - How relevant is the target location to the disease biology?
    - How it the target expression altered in human disease?
    - How is the target involved in the physiological process relevant to the disease?
    - Which phenotypes and genotypes were identified for the target?
    - How is the genetic link between the target and the disease?
    - Describe the evidence provided in clinics or by tools acting on the pathway where the target is involved.
    - Which kind of target modulation is required to treat the disease?

    - Challenges for the drug discovery program related to the target.
    - Check the following idea for details on small molecule compounds: {idea}.
    - Is a 'information driven approach' (IDA) strategy based on available small molecules possible?
        - Which small molecular modulators of the target known?
        - Which inhibitors, antagonists, agonists, negative allosteric modulators (NAM), positive allosteric modulators (PAM) are required for target modulation in the given disease? 
    - Which patients would respond the therapy?
    - Is the proposed mode of action on the target desirable and commercially viable in a clinical setting?
    - What are advantages and disadvantages of different therapeutic modalities (antibodies, small molecules, antisense oligonucleotides, PROTACs, molecular glue, peptide macrocycles, and so on) for tackling the target?

    - Alternative indications:
    - Describe alternative indication for modulators of the target and explain why.

    Task 2: Develop a target assessment strategy for {target} in {disease} in maximal 500 words.

    Outline a 1-year Target Assessment (TA) to Lead Identification (LI) plan. Describe High Level TA-LI plans.
    - Make an emphasis on key inflection points that will inform the feasibility of the project. 
    - Address status of in-vitro platforms, translational in vivo models (mechanistic models, not necessarily so called 'disease models')
    and describe what needs to be established. Elaborate on tractability and major challenges for advancement in a drug discovery portfolio.
    - Discuss potential biomarkers and readouts for efficacy and target engagement.

    Task 3: Safety assessment
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

    Provide the corresponding literature references in the format (First author, Journal, Year, Volume, Issue, Pages, DOI). If any are missing or not available, please mark as N/A.

    Let's work this out in a step by step way to be sure we have the right answer."""
    
    return prompt

target = "Gamma secretase"
disease = "Alzheimer's disease"
mode_of_action = """Gamma secretase is a multi-subunit protease complex that cleaves type I transmembrane proteins, including the amyloid precursor protein (APP) leading to the generation of amyloid-beta (Aβ) peptides."""
context =""" Aβ is a family of secreted peptides generated from the sequential cleavages of the type 1 membrane protein APP by beta-secretase (BACE) and gamma-secretase (GSEC), respectively. BACE cleaves APP in the luminal domain, releasing the N-terminal soluble APPβ domain and leaving the C-terminal fragment, APP-CTF, which remains in the membrane. Subsequently, the APP-CTF is recruited to GSEC, a complex comprising four subunits, including PS, which harbors the active site. GSEC first cuts APP-CTF at the epsilon-cleavage site located close to the inner leaflet of the membrane. This cleavage event produces either Aβ48 or Aβ49 and the APP intracellular domain (AICD). The membrane-retained Aβ48 or Aβ49 is then further processed by GSEC in a continuous cascade of proteolytical events at every third of fourth amino acid, where the N-terminal product of each reaction becomes the substrate for the next GSEC cleavage event. Accordingly, GSEC processes APP-CTF along two main product lines, Aβ49 → 46 → 43 → 40 → 37… and Aβ48 → 45 → 42 → 38…, respectively (Takami et al., 2009; Matsumura et al., 2014; Olsson et al., 2014). During this processing cascade, Aβ43 and shorter Aβ peptides stochastically escape further processing by GSEC and are released into the extracellular space. As a result, Aβ peptides varying from 30 to 43 amino acids in length are secreted into the extracellular space. Among all secreted Aβ, Aβ40 is the most abundant in human CSF, followed by Aβ38, Aβ42, and Aβ37 (Liu et al., 2022). In cognitively normal individuals, Aβ42 and Aβ43 represent a smaller portion of the total secreted Aβ (Liu et al., 2022). These longer forms of Aβ seed the formation of Aβ-amyloid aggregates, a key step in the formation of amyloid plaques (Veugelen et al., 2016), as illustrated in Figure 1. Aβ42, which is produced in higher amounts than Aβ43, is the most abundant Aβ in amyloid plaques (Welander et al., 2009).
""" 

idea = "Develop an inhibitor of the neurokinin 3 receptor (NK3R) GPCR for the treatment of vasomotor symptoms due to menopause."

new_prompt = prompt_to_use(target=target, disease=disease, mode_of_action=mode_of_action, context=context, idea=idea)
print(new_prompt)