from __future__ import annotations

import ast
from typing import Any

import pandas as pd
import plotly.express as px
import ipywidgets as widgets
from IPython.display import display


def _normalize_pgs_ids(value: Any) -> list[str]:
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
    out = df.copy()

    def split_categories(x):
        if x is None:
            return []
        if isinstance(x, float) and pd.isna(x):
            return []
        if isinstance(x, list):
            vals = []
            for item in x:
                if item is None:
                    continue
                item = str(item).strip()
                if not item:
                    continue
                vals.extend([p.strip() for p in item.split(sep) if p.strip()])
            return vals
        if isinstance(x, str):
            return [p.strip() for p in x.split(sep) if p.strip()]
        return [str(x).strip()]

    out[trait_col] = out[trait_col].apply(split_categories)
    out = out.explode(trait_col, ignore_index=True)
    out = out[out[trait_col].notna() & (out[trait_col] != "")].copy()
    return out


def trait_category_pie_with_filter(
    df: pd.DataFrame,
    trait_col: str = "trait_categories",
    label_col: str = "label",
    pgs_col: str = "associated_pgs_ids",
    sep: str = "; ",
    max_rows: int = 500,
):
    df_plot = preprocess_trait_categories(df, trait_col=trait_col, sep=sep).copy()

    if label_col not in df_plot.columns:
        raise ValueError(f"Missing column: {label_col}")
    if pgs_col not in df_plot.columns:
        raise ValueError(f"Missing column: {pgs_col}")

    df_plot[pgs_col] = df_plot[pgs_col].apply(_normalize_pgs_ids)
    df_plot["associated_pgs_ids_str"] = df_plot[pgs_col].apply(lambda x: ", ".join(x))

    counts = (
        df_plot[trait_col]
        .value_counts()
        .rename_axis(trait_col)
        .reset_index(name="count")
        .sort_values("count", ascending=False)
        .reset_index(drop=True)
    )

    fig = px.pie(
        counts,
        names=trait_col,
        values="count",
        title="Traits by category",
    )
    fig.update_traces(textinfo="label+percent")
    fig.update_layout(height=700, margin=dict(l=20, r=20, t=60, b=20))

    category_options = ["-- select category --"] + counts[trait_col].tolist()
    category_dd = widgets.Dropdown(
        options=category_options,
        value="-- select category --",
        description="Category:",
        layout=widgets.Layout(width="100%"),
    )

    summary_html = widgets.HTML("<i>Select a category from the dropdown.</i>")
    table_out = widgets.Output(
        layout=widgets.Layout(
            border="1px solid #ddd",
            height="620px",
            overflow_y="auto",
            padding="6px",
        )
    )

    def render_table(category: str | None):
        with table_out:
            table_out.clear_output()

            if not category or category == "-- select category --":
                summary_html.value = "<i>Select a category from the dropdown.</i>"
                display(pd.DataFrame(columns=[label_col, pgs_col]))
                return

            subset = (
                df_plot.loc[df_plot[trait_col] == category, [label_col, "associated_pgs_ids_str"]]
                .rename(columns={"associated_pgs_ids_str": pgs_col})
                .drop_duplicates()
                .sort_values(label_col)
                .reset_index(drop=True)
            )

            n = len(subset)
            summary_html.value = f"<b>Category:</b> {category}<br><b>Traits:</b> {n}"

            if n > max_rows:
                display(subset.head(max_rows))
                print(f"Showing first {max_rows} of {n} rows.")
            else:
                display(subset)

    def _on_change(change):
        if change["name"] == "value":
            render_table(change["new"])

    category_dd.observe(_on_change, names="value")

    left = widgets.VBox(
        [widgets.HTML(fig.to_html(include_plotlyjs="cdn"))],
        layout=widgets.Layout(width="55%"),
    )

    right = widgets.VBox(
        [widgets.HTML("<h3>Selected category</h3>"), category_dd, summary_html, table_out],
        layout=widgets.Layout(width="45%"),
    )

    app = widgets.HBox([left, right], layout=widgets.Layout(width="100%"))
    render_table(None)
    return app