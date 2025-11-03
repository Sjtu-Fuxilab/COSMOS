# Detailed Methodology

## Cox Proportional Hazards

**Model:** `h(t|X) = h₀(t) × exp(β × Phase)`

**Patient Aggregation:** `patient_embedding = mean(slide_embeddings)`

**Tests:** Cox LRT, Wald test, Log-rank, C-index

**Power:** 1,362 events (3.1× required for 90% power)

---

## Feature Extraction

1. Extract 224×224 patches from WSI
2. ResNet50 → 2,048-D features
3. Slide-level: mean pooling
4. Patient-level: mean across slides

---

## Implementation

- lifelines v0.27.4
- PyTorch v1.12.0  
- OpenSlide v1.2.0
- Fixed seeds: `np.random.seed(42)`
