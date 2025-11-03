# Installation Guide

## Quick Install

```bash
pip install -r requirements.txt
```

---

## Platform-Specific

### Linux
```bash
sudo apt-get install openslide-tools python3-openslide
pip install -r requirements.txt
```

### macOS
```bash
brew install openslide
pip install -r requirements.txt
```

### Windows
1. Download OpenSlide from: https://openslide.org/download/
2. Add to PATH
3. `pip install -r requirements.txt`

---

## Troubleshooting

### CUDA Out of Memory
```bash
python full_validation.py --batch_size 64
```

### Verify Installation
```python
import torch, openslide, lifelines
print("✅ All packages installed")
```
