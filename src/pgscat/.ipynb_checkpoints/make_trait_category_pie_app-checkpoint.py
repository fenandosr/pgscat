from __future__ import annotations

import ast
from typing import Any

import pandas as pd
import plotly.graph_objects as go
import ipywidgets as widgets
from IPython.display import display


def _normalize_pgs_ids(value: Any) -> list[str]:
    """
    Normaliza associated_pgs_ids a list[str].
    Soporta:
    - list
    - tuple
    - None / NaN
    - string tipo "['PGS0001', 'PGS0002']"
    - string simple
    """
    if value is None:
        return []

    if isinstance(value, float) and pd.isna(value):
        return []

    if isinstance(value, (list, tuple)):
        return [str(x) for x in value]

    if isinstance(value, str):
        s = value.strip()
        if not s:
            return []

        # Intenta parsear listas serializadas
        if s.startswith("[") and s.endswith("]"):
            try:
                parsed = ast.literal_eval(s)
                if isinstance(parsed, (list, tuple)):
                    return [str(x) for x in parsed]
            except Exception:
                pass

        return [s]

    return [str(value)]


def preprocess_trait_categories(
    df: pd.DataFrame,
    trait_col: str = "trait_categories",
    sep: str = "; ",
) -> pd.DataFrame:
    """
    Separa categorías múltiples en trait_categories y duplica registros.

    Ejemplo:
      'Cardiovascular disease; Metabolic disease'
    ->
      dos filas, una por categoría.
    """
    out = df.copy()

    # Asegurar string y manejar missing
    out[trait_col] = out[trait_col].fillna("").astype(str).str.strip()

    # split
    out[trait_col] = out[trait_col].apply(
        lambda x: [p.strip() for p in x.split(sep) if p.strip()] if x else []
    )

    # explode
    out = out.explode(trait_col, ignore_index=True)

    # eliminar vacíos por si hubo rows sin categoría
    out = out[out[trait_col].notna() & (out[trait_col] != "")].copy()

    return out


def make_trait_category_pie_app(
    df: pd.DataFrame,
    trait_col: str = "trait_categories",
    label_col: str = "label",
    pgs_col: str = "associated_pgs_ids",
    sep: str = "; ",
    max_rows_in_panel: int = 500,
) -> widgets.Widget:
    """
    Crea una app interactiva:
    - pie chart de trait_categories
    - click en una categoría -> panel derecho con label y associated_pgs_ids
    """
    df_plot = preprocess_trait_categories(df, trait_col=trait_col, sep=sep).copy()

    if label_col not in df_plot.columns:
        raise ValueError(f"Missing required column: {label_col}")

    if pgs_col not in df_plot.columns:
        raise ValueError(f"Missing required column: {pgs_col}")

    df_plot[pgs_col] = df_plot[pgs_col].apply(_normalize_pgs_ids)
    df_plot["associated_pgs_ids_str"] = df_plot[pgs_col].apply(lambda x: ", ".join(x))

    # Conteo por categoría
    counts = (
        df_plot[trait_col]
        .value_counts(dropna=False)
        .rename_axis(trait_col)
        .reset_index(name="count")
        .sort_values("count", ascending=False)
        .reset_index(drop=True)
    )

    # Pie chart
    fig = go.FigureWidget(
        data=[
            go.Pie(
                labels=counts[trait_col],
                values=counts["count"],
                sort=False,
                direction="clockwise",
                textinfo="label+percent",
                hovertemplate="<b>%{label}</b><br>count=%{value}<extra></extra>",
            )
        ]
    )

    fig.update_layout(
        title="Traits by category",
        height=700,
        margin=dict(l=20, r=20, t=60, b=20),
    )

    # Widgets del panel derecho
    title_html = widgets.HTML("<h3>Selected category</h3>")
    summary_html = widgets.HTML("<i>Click a pie slice to see matching traits.</i>")

    output_table = widgets.Output(
        layout=widgets.Layout(
            border="1px solid #ddd",
            height="620px",
            overflow_y="auto",
            padding="6px",
        )
    )

    reset_button = widgets.Button(
        description="Reset selection",
        button_style="",
        tooltip="Clear selection and show instructions",
        icon="refresh",
    )

    def _render_category(category: str | None) -> None:
        with output_table:
            output_table.clear_output()

            if category is None:
                summary_html.value = "<i>Click a pie slice to see matching traits.</i>"
                display(pd.DataFrame(columns=[label_col, pgs_col]))
                return

            subset = (
                df_plot.loc[df_plot[trait_col] == category, [label_col, "associated_pgs_ids_str"]]
                .rename(columns={"associated_pgs_ids_str": pgs_col})
                .drop_duplicates()
                .reset_index(drop=True)
            )

            n_rows = len(subset)
            summary_html.value = (
                f"<b>Category:</b> {category}<br>"
                f"<b>Traits:</b> {n_rows}"
            )

            if n_rows > max_rows_in_panel:
                display(subset.head(max_rows_in_panel))
                display(
                    widgets.HTML(
                        f"<p><i>Showing first {max_rows_in_panel} of {n_rows} rows.</i></p>"
                    )
                )
            else:
                display(subset)

    def _handle_click(trace, points, state):
        if not points.point_inds:
            return

        idx = points.point_inds[0]
        category = counts.iloc[idx][trait_col]
        _render_category(category)

    def _reset_clicked(_):
        _render_category(None)

    fig.data[0].on_click(_handle_click)
    reset_button.on_click(_reset_clicked)

    right_panel = widgets.VBox(
        [title_html, summary_html, reset_button, output_table],
        layout=widgets.Layout(width="45%"),
    )

    left_panel = widgets.VBox(
        [fig],
        layout=widgets.Layout(width="55%"),
    )

    app = widgets.HBox(
        [left_panel, right_panel],
        layout=widgets.Layout(width="100%"),
    )

    _render_category(None)
    return app