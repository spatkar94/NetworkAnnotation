# Network annotation
A chief goal of systems biology is the reconstruction of large-scale executable models of
cellular processes of interest. While accurate continuous models are still beyond reach, a powerful
alternative is to learn a logical model of the processes under study, which predicts the logical state of
any node of the model as a Boolean function of its incoming nodes. Key to learning such models is the
functional annotation of the underlying physical interactions (protein-protein or protein-DNA) with activation/repression (sign) effects. Such annotations are pretty common for a few well studied biological pathways. Here, we developed a novel optimization framework that uses different plausible models
of cellular signaling to predict the signs of yet unannotated physical interactions from large-scale gene-knockout experiments on a system-wide level. 

Currently, the framework is implemented on physical interactions of budding yeast. (saccharomyces cerevisiae)

# Requirements
- Gurobi (>=7.0.2) (Click [here](https://www.gurobi.com/documentation/7.5/quickstart_linux/the_gurobi_python_interfac.html) for installation instructions).

- Python (>=2.7.1).

Run the following from the command line to install remaining pre-requisite packages.
```
    cd /my_path/NetworkAnnotation/
    pip install -r requirements.txt
```
# Usage
Run the following from the command line
```
    cd /my_path/NetworkAnnotation/
    python assign.py -v AllSP -i reimand.txt -r 10 -f 5

```

-v: The signaling model used.
- AP (stands for A path. https://dl.acm.org/citation.cfm?id=2415282)
- ASP (stands for A Shortest Path) 
- AdirSP (stands for A directed Shortest Path)
- AllSP (stands for All Shortest Paths)

-i: input file consisting of knock-out pairs inferred from gene knockout experiments via differential expression analysis. (See reimand.txt)

-r: number of repeated optimizations to measure sign confidence scores. A score near 0 implies a positive sign prediction with high confidence, a score near 1 implies a negative sign prediction with high confidence. A score near 0.5 either implies ambiguity of sign or a lack of a sign. (Note, not all physical interactions are expected to possess a sign)

-f: number of folds of cross-validation

# Output
A table consisting of the actual signs and predicted sign confidence scores of all pre-annotated physical interactions along with their type (kpi: kinase-substrate/phosphatase-substrate interactions, pdi: protein DNA/regulatory interactions, KEGG: experimentally annotated physical interactions from the KEGG database, http://www.genome.jp/kegg/pathway.html). See output.txt.

# Link to publication 
Currently in press.
