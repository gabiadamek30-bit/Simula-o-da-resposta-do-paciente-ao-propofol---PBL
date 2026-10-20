"""
NeoPKPD Visualization Themes

Professional color palettes and styling for PK/PD visualizations.
"""

from __future__ import annotations

from typing import Dict, List, Optional

# Global theme setting
_CURRENT_THEME = "neopkpd"


# ============================================================================
# Color Palettes
# ============================================================================

# NeoPKPD brand colors - professional, accessible palette
NEOPKPD_COLORS = {
    "primary": "#2E86AB",      # Blue - main color
    "secondary": "#A23B72",    # Magenta - secondary
    "accent": "#F18F01",       # Orange - highlights
    "success": "#4CAF50",      # Green - success/reference
    "warning": "#FFC107",      # Yellow - warnings
    "error": "#E94560",        # Red - errors/high values
    "info": "#17A2B8",         # Cyan - info
    "dark": "#2C3E50",         # Dark gray - text
    "light": "#ECF0F1",        # Light gray - backgrounds
    "white": "#FFFFFF",        # White
}

# Color sequences for multi-series plots
COLOR_SEQUENCES = {
    "neopkpd": [
        "#2E86AB",  # Blue
        "#E94560",  # Red
        "#4CAF50",  # Green
        "#F18F01",  # Orange
        "#A23B72",  # Magenta
        "#17A2B8",  # Cyan
        "#9C27B0",  # Purple
        "#795548",  # Brown
    ],
    "clinical": [
        "#1f77b4",  # Blue
        "#ff7f0e",  # Orange
        "#2ca02c",  # Green
        "#d62728",  # Red
        "#9467bd",  # Purple
        "#8c564b",  # Brown
        "#e377c2",  # Pink
        "#7f7f7f",  # Gray
    ],
    "grayscale": [
        "#000000",  # Black
        "#333333",  # Dark gray
        "#666666",  # Medium gray
        "#999999",  # Light gray
        "#CCCCCC",  # Very light gray
    ],
}

# Confidence interval colors (transparent versions)
CI_COLORS = {
    "90": "rgba(46, 134, 171, 0.2)",   # 90% CI - light blue
    "50": "rgba(46, 134, 171, 0.4)",   # 50% CI - medium blue
    "median": "#2E86AB",               # Median line
}


# ============================================================================
# Theme Definitions
# ============================================================================

THEMES = {
    "neopkpd": {
        "colors": COLOR_SEQUENCES["neopkpd"],
        "background": "#FFFFFF",
        "text_color": "#2C3E50",
        "grid_color": "#E0E0E0",
        "font_family": "sans-serif",
        "font_size": 12,
        "title_size": 14,
        "line_width": 1.5,
        "marker_size": 6,
        "alpha_main": 1.0,
        "alpha_ribbon": 0.3,
        "alpha_spaghetti": 0.4,
    },
    "clinical": {
        "colors": COLOR_SEQUENCES["clinical"],
        "background": "#FFFFFF",
        "text_color": "#000000",
        "grid_color": "#CCCCCC",
        "font_family": "Arial",
        "font_size": 11,
        "title_size": 13,
        "line_width": 1.0,
        "marker_size": 5,
        "alpha_main": 1.0,
        "alpha_ribbon": 0.25,
        "alpha_spaghetti": 0.3,
    },
    "presentation": {
        "colors": COLOR_SEQUENCES["neopkpd"],
        "background": "#FFFFFF",
        "text_color": "#2C3E50",
        "grid_color": "#E0E0E0",
        "font_family": "Arial",
        "font_size": 14,
        "title_size": 18,
        "line_width": 2.5,
        "marker_size": 8,
        "alpha_main": 1.0,
        "alpha_ribbon": 0.3,
        "alpha_spaghetti": 0.5,
    },
    "publication": {
        "colors": COLOR_SEQUENCES["grayscale"],
        "background": "#FFFFFF",
        "text_color": "#000000",
        "grid_color": "#CCCCCC",
        "font_family": "Times New Roman",
        "font_size": 10,
        "title_size": 12,
        "line_width": 1.0,
        "marker_size": 4,
        "alpha_main": 1.0,
        "alpha_ribbon": 0.2,
        "alpha_spaghetti": 0.3,
    },
}


# ============================================================================
# Theme Functions
# ============================================================================

def available_themes() -> List[str]:
    """
    Get list of available themes.

    Returns:
        List of theme names
    """
    return list(THEMES.keys())


def get_theme() -> str:
    """
    Get the current theme name.

    Returns:
        Current theme name
    """
    return _CURRENT_THEME


def set_theme(theme: str) -> None:
    """
    Set the visualization theme.

    Args:
        theme: Theme name

    Raises:
        ValueError: If theme is not available

    Example:
        >>> from neopkpd import viz
        >>> viz.set_theme("presentation")
    """
    global _CURRENT_THEME

    if theme not in THEMES:
        raise ValueError(
            f"Theme '{theme}' is not available. "
            f"Available themes: {list(THEMES.keys())}"
        )

    _CURRENT_THEME = theme


def get_theme_config() -> Dict:
    """
    Get the current theme configuration.

    Returns:
        Theme configuration dictionary
    """
    return THEMES[_CURRENT_THEME].copy()


def get_color(index: int) -> str:
    """
    Get color by index from current theme.

    Args:
        index: Color index (wraps around)

    Returns:
        Color hex string
    """
    colors = THEMES[_CURRENT_THEME]["colors"]
    return colors[index % len(colors)]


def get_colors(n: int) -> List[str]:
    """
    Get n colors from current theme.

    Args:
        n: Number of colors needed

    Returns:
        List of color hex strings
    """
    colors = THEMES[_CURRENT_THEME]["colors"]
    return [colors[i % len(colors)] for i in range(n)]
