import pytest
from pathlib import Path

from pgscat.common_utils import file_delim_meta_header_data


def write_file(tmp_path: Path, content: str) -> Path:
    p = tmp_path / "test.txt"
    p.write_text(content)
    return p


def test_basic_commented_file(tmp_path):
    content = """\
## comment
#weight_type=beta
#source=test

chr_name\tchr_position\teffect_allele\tother_allele\teffect_weight
1\t12345\tA\tG\t0.1
1\t12346\tC\tT\t0.2
"""
    path = write_file(tmp_path, content)

    delim, meta, header, data = file_delim_meta_header_data(
        str(path),
        meta_line_startswith="#",
        meta_kv_separator="=",
    )

    assert delim == "\t"
    assert header == [
        "chr_name",
        "chr_position",
        "effect_allele",
        "other_allele",
        "effect_weight",
    ]
    assert meta["line1"] == "comment"
    assert meta["weight_type"] == "beta"
    assert meta["source"] == "test"
    assert data[0] == ["1", "12345", "A", "G", "0.1"]


def test_csv_delimiter(tmp_path):
    content = """\
## comment
#weight_type=beta
chr,pos,ea,oa,w
1,100,A,G,0.5
"""
    path = write_file(tmp_path, content)

    delim, meta, header, data = file_delim_meta_header_data(
        str(path),
        meta_line_startswith="#",
        skip_meta_line_startswith="##",
        meta_kv_separator="=",
    )
    assert delim == ","
    assert header == ["chr", "pos", "ea", "oa", "w"]
    assert meta["weight_type"] == "beta"
    assert len(meta) == 1
    assert data[0] == ["1", "100", "A", "G", "0.5"]


def test_forced_delimiter(tmp_path):
    content = """\
#k=v
a|b|c
1|2|3
"""
    path = write_file(tmp_path, content)

    delim, meta, header, data = file_delim_meta_header_data(
        str(path),
        forced_delimiter="|",
        skip_meta_line_startswith="#",
    )

    assert delim == "|"
    assert len(meta) == 0
    assert header == ["a", "b", "c"]
    assert data == [["1", "2", "3"]]


def test_no_meta(tmp_path):
    content = """\
a\tb\tc
1\t2\t3
"""
    path = write_file(tmp_path, content)

    delim, meta, header, data = file_delim_meta_header_data(str(path))
    assert delim == "\t"
    assert meta == {}
    assert header == ["a", "b", "c"]
    assert data == [["1", "2", "3"]]


def test_empty_lines_and_comments(tmp_path):
    content = """\

## ignore
#foo=bar

a\tb\tc
1\t2\t3

"""
    path = write_file(tmp_path, content)

    delim, meta, header, data = file_delim_meta_header_data(
        str(path),
        meta_line_startswith="#",
        skip_meta_line_startswith="##",
        meta_kv_separator="=",
    )
    assert delim == "\t"
    assert header == ["a", "b", "c"]
    assert meta["foo"] == "bar"
    assert len(meta) == 1
    assert data == [["1", "2", "3"]]


def test_no_data_raises(tmp_path):
    content = """\
## only comments
#k=v
"""
    path = write_file(tmp_path, content)

    with pytest.raises(ValueError):
        file_delim_meta_header_data(
            str(path),
            meta_line_startswith="#",
        )


def test_skip_data(tmp_path):
    content = """\
## comment
#weight_type=beta
chr,pos,ea,oa,w
1,100,A,G,0.5
"""
    path = write_file(tmp_path, content)

    delim, meta, header, data = file_delim_meta_header_data(
        str(path),
        meta_line_startswith="#",
        skip_meta_line_startswith="##",
        meta_kv_separator="=",
        skip_data=True,
    )
    assert data == []
