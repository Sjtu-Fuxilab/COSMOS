# Detailed Methodology

## Cox Proportional Hazards Regression

### Model Specification

```
h(t|X) = h₀(t) × exp(β × Phase)
```

Where:
- h(t|X) = hazard function
- Phase = binary (1 = early ≤24mo, 0 = late >24mo)
- HR = exp(β)

### Patient-Level Aggregation

**Problem:** Multiple slides per patient violate independence.

**Solution:**
```python
patient_embedding = mean(slide_embeddings)
```

### Statistical Tests

1. Cox Likelihood Ratio Test
2. Wald Test (95% CI)
3. Log-Rank Test
4. Harrell's C-index

### Sample Size

Required: ~440 events (90% power for HR=1.5)  
Actual: **1,362 events** (3.1× requirement)

---

## Feature Extraction

1. **Patch Extraction:** 224×224 from WSIs
2. **ResNet50:** 2,048-D features per patch
3. **Slide Aggregation:** Mean pooling
4. **Patient Aggregation:** Mean across slides

---

## Implementation

**Libraries:**
- lifelines v0.27.4
- PyTorch v1.12.0
- OpenSlide v1.2.0

**Reproducibility:**
```python
np.random.seed(42)
torch.manual_seed(42)
```
