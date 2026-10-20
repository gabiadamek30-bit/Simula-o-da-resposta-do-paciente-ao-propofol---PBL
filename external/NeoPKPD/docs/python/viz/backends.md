# Backends & Themes

Configure visualization backends and styling for NeoPKPD plots.

---

## Overview

NeoPKPD supports dual visualization backends:

- **Matplotlib**: Static, publication-quality figures
- **Plotly**: Interactive, web-embeddable plots

```python
from neopkpd import viz

# Set backend globally
viz.set_backend("matplotlib")  # or "plotly"

# Check current backend
print(viz.get_backend())

# List available backends
print(viz.available_backends())  # ["matplotlib", "plotly"]
```

---

## Backend Selection

### Matplotlib (Default)

Best for publication-quality static figures:

```python
viz.set_backend("matplotlib")

fig = viz.plot_conc_time(result)
fig.savefig("plot.pdf", bbox_inches="tight")  # Vector format
fig.savefig("plot.png", dpi=300)              # Raster format
```

**Advantages:**
- Publication-quality output
- Vector format support (PDF, SVG, EPS)
- Fine-grained customization
- Familiar API for scientists

### Plotly

Best for interactive exploration:

```python
viz.set_backend("plotly")

fig = viz.plot_conc_time(result)
fig.write_html("plot.html")           # Interactive HTML
fig.write_image("plot.png", scale=2)  # Static image (requires kaleido)
```

**Advantages:**
- Interactive zoom/pan
- Hover tooltips
- Web embedding
- Animation support

---

## Themes

### Available Themes

```python
# Set theme
viz.set_theme("neopkpd")      # Default professional theme
viz.set_theme("publication")    # Minimal for publications
viz.set_theme("presentation")   # Bold for slides

# List available themes
print(viz.available_themes())
```

### Theme Properties

| Theme | Use Case | Font Size | Line Width |
|-------|----------|-----------|------------|
| `neopkpd` | General use | Medium | Medium |
| `publication` | Journal figures | Small | Thin |
| `presentation` | Slides | Large | Thick |

---

## Color Palette

```python
# Access color palette
colors = viz.NEOPKPD_COLORS

print(colors)
# {
#     "primary": "#3498DB",
#     "secondary": "#2ECC71",
#     "accent": "#E74C3C",
#     "neutral": "#95A5A6",
#     "dark": "#2C3E50",
#     "light": "#ECF0F1"
# }

# Use in custom plots
import matplotlib.pyplot as plt
plt.plot(x, y, color=colors["primary"])
```

---

## Per-Function Backend Override

Override backend for individual function calls:

```python
# Global backend is matplotlib
viz.set_backend("matplotlib")

# But use plotly for this specific plot
fig = viz.plot_conc_time(result, backend="plotly")
```

---

## Saving Figures

### Matplotlib

```python
viz.set_backend("matplotlib")
fig = viz.plot_conc_time(result)

# PNG (raster)
fig.savefig("plot.png", dpi=300, bbox_inches="tight")

# PDF (vector)
fig.savefig("plot.pdf", bbox_inches="tight")

# SVG (vector, web-friendly)
fig.savefig("plot.svg", bbox_inches="tight")

# Using save_path parameter
fig = viz.plot_conc_time(result, save_path="plot.png")
```

### Plotly

```python
viz.set_backend("plotly")
fig = viz.plot_conc_time(result)

# Interactive HTML
fig.write_html("plot.html")

# Static image (requires kaleido: pip install kaleido)
fig.write_image("plot.png", scale=2)
fig.write_image("plot.pdf")
fig.write_image("plot.svg")
```

---

## Custom Styling

### Matplotlib Customization

```python
import matplotlib.pyplot as plt

# Custom rcParams
plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 12,
    'axes.linewidth': 1.5,
    'axes.labelsize': 14,
    'xtick.labelsize': 11,
    'ytick.labelsize': 11,
    'legend.fontsize': 11,
})

# Apply to plot
fig = viz.plot_conc_time(result)
```

### Plotly Customization

```python
fig = viz.plot_conc_time(result, backend="plotly")

# Update layout
fig.update_layout(
    font=dict(family="Arial", size=14),
    plot_bgcolor="white",
    paper_bgcolor="white"
)

# Update traces
fig.update_traces(line=dict(width=2))
```

---

## See Also

- [PK Plots](pk.md) - Concentration-time visualization
- [VPC Plots](vpc.md) - Visual predictive checks
- [Visualization Index](index.md) - All visualization functions
