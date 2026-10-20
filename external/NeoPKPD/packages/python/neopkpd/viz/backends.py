"""
NeoPKPD Visualization Backends

Provides infrastructure for multiple plotting backends (matplotlib, plotly).
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from enum import Enum
from typing import Any, Dict, List, Optional, Tuple, Union

# Global backend setting
_CURRENT_BACKEND = "matplotlib"


class PlotBackend(str, Enum):
    """Available plotting backends."""
    MATPLOTLIB = "matplotlib"
    PLOTLY = "plotly"


def available_backends() -> List[str]:
    """
    Get list of available plotting backends.

    Returns:
        List of backend names
    """
    available = []

    try:
        import matplotlib
        available.append("matplotlib")
    except ImportError:
        pass

    try:
        import plotly
        available.append("plotly")
    except ImportError:
        pass

    return available


def get_backend() -> str:
    """
    Get the current plotting backend.

    Returns:
        Current backend name
    """
    return _CURRENT_BACKEND


def set_backend(backend: str) -> None:
    """
    Set the plotting backend.

    Args:
        backend: Backend name ("matplotlib" or "plotly")

    Raises:
        ValueError: If backend is not available

    Example:
        >>> from neopkpd import viz
        >>> viz.set_backend("plotly")
    """
    global _CURRENT_BACKEND

    backend = backend.lower()
    available = available_backends()

    if backend not in available:
        raise ValueError(
            f"Backend '{backend}' is not available. "
            f"Available backends: {available}. "
            f"Install with: pip install {backend}"
        )

    _CURRENT_BACKEND = backend


def _get_plotter():
    """Get the current backend plotter."""
    backend = get_backend()

    if backend == "matplotlib":
        return MatplotlibPlotter()
    elif backend == "plotly":
        return PlotlyPlotter()
    else:
        raise ValueError(f"Unknown backend: {backend}")


# ============================================================================
# Backend Abstract Interface
# ============================================================================

class BasePlotter(ABC):
    """Abstract base class for plotting backends."""

    @abstractmethod
    def create_figure(
        self,
        figsize: Tuple[float, float] = (10, 6),
        title: Optional[str] = None,
    ) -> Any:
        """Create a new figure."""
        pass

    @abstractmethod
    def line_plot(
        self,
        fig: Any,
        x: List[float],
        y: List[float],
        label: Optional[str] = None,
        color: Optional[str] = None,
        linestyle: str = "-",
        linewidth: float = 1.5,
        alpha: float = 1.0,
    ) -> Any:
        """Add a line to the figure."""
        pass

    @abstractmethod
    def scatter_plot(
        self,
        fig: Any,
        x: List[float],
        y: List[float],
        label: Optional[str] = None,
        color: Optional[str] = None,
        marker: str = "o",
        size: float = 30,
        alpha: float = 1.0,
    ) -> Any:
        """Add scatter points to the figure."""
        pass

    @abstractmethod
    def fill_between(
        self,
        fig: Any,
        x: List[float],
        y_lower: List[float],
        y_upper: List[float],
        color: Optional[str] = None,
        alpha: float = 0.3,
        label: Optional[str] = None,
    ) -> Any:
        """Add filled region between two lines."""
        pass

    @abstractmethod
    def set_labels(
        self,
        fig: Any,
        xlabel: Optional[str] = None,
        ylabel: Optional[str] = None,
        title: Optional[str] = None,
    ) -> None:
        """Set axis labels and title."""
        pass

    @abstractmethod
    def set_log_scale(
        self,
        fig: Any,
        x: bool = False,
        y: bool = False,
    ) -> None:
        """Set logarithmic scale for axes."""
        pass

    @abstractmethod
    def add_legend(
        self,
        fig: Any,
        location: str = "best",
    ) -> None:
        """Add legend to figure."""
        pass

    @abstractmethod
    def finalize(self, fig: Any) -> Any:
        """Finalize and return the figure."""
        pass


# ============================================================================
# Matplotlib Backend
# ============================================================================

class MatplotlibPlotter(BasePlotter):
    """Matplotlib backend implementation."""

    def __init__(self):
        import matplotlib.pyplot as plt
        self.plt = plt

    def create_figure(
        self,
        figsize: Tuple[float, float] = (10, 6),
        title: Optional[str] = None,
    ) -> Any:
        fig, ax = self.plt.subplots(figsize=figsize)
        if title:
            ax.set_title(title)
        return {"fig": fig, "ax": ax}

    def line_plot(
        self,
        fig: Any,
        x: List[float],
        y: List[float],
        label: Optional[str] = None,
        color: Optional[str] = None,
        linestyle: str = "-",
        linewidth: float = 1.5,
        alpha: float = 1.0,
    ) -> Any:
        ax = fig["ax"]
        ax.plot(x, y, label=label, color=color, linestyle=linestyle,
                linewidth=linewidth, alpha=alpha)
        return fig

    def scatter_plot(
        self,
        fig: Any,
        x: List[float],
        y: List[float],
        label: Optional[str] = None,
        color: Optional[str] = None,
        marker: str = "o",
        size: float = 30,
        alpha: float = 1.0,
    ) -> Any:
        ax = fig["ax"]
        ax.scatter(x, y, label=label, color=color, marker=marker,
                   s=size, alpha=alpha)
        return fig

    def fill_between(
        self,
        fig: Any,
        x: List[float],
        y_lower: List[float],
        y_upper: List[float],
        color: Optional[str] = None,
        alpha: float = 0.3,
        label: Optional[str] = None,
    ) -> Any:
        ax = fig["ax"]
        ax.fill_between(x, y_lower, y_upper, color=color, alpha=alpha, label=label)
        return fig

    def set_labels(
        self,
        fig: Any,
        xlabel: Optional[str] = None,
        ylabel: Optional[str] = None,
        title: Optional[str] = None,
    ) -> None:
        ax = fig["ax"]
        if xlabel:
            ax.set_xlabel(xlabel)
        if ylabel:
            ax.set_ylabel(ylabel)
        if title:
            ax.set_title(title)

    def set_log_scale(
        self,
        fig: Any,
        x: bool = False,
        y: bool = False,
    ) -> None:
        ax = fig["ax"]
        if x:
            ax.set_xscale("log")
        if y:
            ax.set_yscale("log")

    def add_legend(
        self,
        fig: Any,
        location: str = "best",
    ) -> None:
        ax = fig["ax"]
        ax.legend(loc=location)

    def finalize(self, fig: Any) -> Any:
        fig["fig"].tight_layout()
        return fig["fig"]


# ============================================================================
# Plotly Backend
# ============================================================================

class PlotlyPlotter(BasePlotter):
    """Plotly backend implementation."""

    def __init__(self):
        import plotly.graph_objects as go
        self.go = go

    def create_figure(
        self,
        figsize: Tuple[float, float] = (10, 6),
        title: Optional[str] = None,
    ) -> Any:
        fig = self.go.Figure()
        if title:
            fig.update_layout(title=title)
        # Convert figsize from inches to pixels (approx 100 dpi)
        fig.update_layout(width=figsize[0] * 100, height=figsize[1] * 100)
        return fig

    def line_plot(
        self,
        fig: Any,
        x: List[float],
        y: List[float],
        label: Optional[str] = None,
        color: Optional[str] = None,
        linestyle: str = "-",
        linewidth: float = 1.5,
        alpha: float = 1.0,
    ) -> Any:
        dash_map = {"-": "solid", "--": "dash", ":": "dot", "-.": "dashdot"}
        dash = dash_map.get(linestyle, "solid")

        fig.add_trace(self.go.Scatter(
            x=x, y=y, mode="lines", name=label,
            line=dict(color=color, dash=dash, width=linewidth),
            opacity=alpha,
        ))
        return fig

    def scatter_plot(
        self,
        fig: Any,
        x: List[float],
        y: List[float],
        label: Optional[str] = None,
        color: Optional[str] = None,
        marker: str = "o",
        size: float = 30,
        alpha: float = 1.0,
    ) -> Any:
        marker_map = {"o": "circle", "s": "square", "^": "triangle-up", "d": "diamond"}
        symbol = marker_map.get(marker, "circle")

        fig.add_trace(self.go.Scatter(
            x=x, y=y, mode="markers", name=label,
            marker=dict(color=color, symbol=symbol, size=size / 5),
            opacity=alpha,
        ))
        return fig

    def fill_between(
        self,
        fig: Any,
        x: List[float],
        y_lower: List[float],
        y_upper: List[float],
        color: Optional[str] = None,
        alpha: float = 0.3,
        label: Optional[str] = None,
    ) -> Any:
        # Upper bound
        fig.add_trace(self.go.Scatter(
            x=x, y=y_upper, mode="lines", line=dict(width=0),
            showlegend=False,
        ))
        # Lower bound with fill
        fig.add_trace(self.go.Scatter(
            x=x, y=y_lower, mode="lines", line=dict(width=0),
            fill="tonexty", fillcolor=color if color else "rgba(0,100,200,0.3)",
            name=label, opacity=alpha,
        ))
        return fig

    def set_labels(
        self,
        fig: Any,
        xlabel: Optional[str] = None,
        ylabel: Optional[str] = None,
        title: Optional[str] = None,
    ) -> None:
        layout_updates = {}
        if xlabel:
            layout_updates["xaxis_title"] = xlabel
        if ylabel:
            layout_updates["yaxis_title"] = ylabel
        if title:
            layout_updates["title"] = title
        fig.update_layout(**layout_updates)

    def set_log_scale(
        self,
        fig: Any,
        x: bool = False,
        y: bool = False,
    ) -> None:
        if x:
            fig.update_xaxes(type="log")
        if y:
            fig.update_yaxes(type="log")

    def add_legend(
        self,
        fig: Any,
        location: str = "best",
    ) -> None:
        fig.update_layout(showlegend=True)

    def finalize(self, fig: Any) -> Any:
        return fig
