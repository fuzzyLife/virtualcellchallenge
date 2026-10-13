---
tags:
- single-cell
- biology
- foundation-model
- in-context-learning
- transcriptomics
- genomics
---

# Stack-Large-Aligned

**Stack** is a large-scale encoder-decoder foundation model for single-cell biology. It introduces a novel tabular attention architecture that enables both intra- and inter-cellular information flow, setting cell-by-gene matrix chunks as the basic input data unit. Through in-context learning, Stack offers substantial performance improvements in generalizing biological effects and enables generation of unseen cell profiles in novel contexts.

## Model Details

| Property | Value |
|----------|-------|
| **Parameters** | 217M |
| **Architecture** | Tabular Attention (alternating cell-wise and gene-wise attention) |
| **Model Size** | Large |
| **Pretraining Data** | Full human scBaseCount (~150M cells) |
| **Pretraining Epochs** | 10 |
| **Alignment Data** | CellxGene 45M + Parse 10M PBMC |
| **Alignment Epochs** | 8 |

## Usage

### Installation
```bash
pip install arc-stack
```

### Download Model
```python
from huggingface_hub import snapshot_download

repo_id = "arcinstitute/Stack-Large-Aligned"
local_dir = "Stack-Large-Aligned"
snapshot_download(repo_id=repo_id, repo_type="model", local_dir=local_dir)
```

For detailed tutorials, see:
- [Prediction Tutorial](https://github.com/ArcInstitute/stack/blob/main/notebooks/tutorial-predict.ipynb)


## Citation

If you use this model, please cite:

- Dong et al., 2026: Stack: In-context modeling of single-cell biology. Preprint. [Paper link](https://www.biorxiv.org/content/10.64898/2026.01.09.698608v1)

## License

This model is released under a non-commercial license.

For model licenses please see MODEL_ACCEPTABLE_USE_POLICY.md, MODEL_LICENSE.md, and LICENSE in [stack Github](https://github.com/ArcInstitute/stack).

## Github Link

[https://github.com/ArcInstitute/stack](https://github.com/ArcInstitute/stack)