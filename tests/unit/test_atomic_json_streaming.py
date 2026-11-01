"""Stream artifact JSON with stable bytes and atomic failure behavior."""
import json
from pathlib import Path

import pytest

from zipper_algebra.cofactor import artifacts


@pytest.mark.parametrize('document', [
    {}, [], {'unicode': 'π𝛼', 'escape': '\n\\"'},
    {'mixed': [1, -2, 1.25, -0.0, True, None, {}, []], 'large': 10**80},
    {'nested': [{'words': list(range(40))} for _ in range(20)]},
])
def test_streaming_writer_preserves_existing_json_bytes(tmp_path, document):
    path = tmp_path/'artifact.json'
    artifacts.write_json_atomic(path, document)
    expected = json.dumps(document, sort_keys=True, indent=2, allow_nan=False)+'\n'
    assert path.read_bytes() == expected.encode('utf-8')
    assert json.loads(path.read_text()) == json.loads(expected)
    assert sorted(p.name for p in tmp_path.iterdir()) == ['artifact.json']


def test_streaming_writer_does_not_build_a_whole_encoded_string(tmp_path, monkeypatch):
    # A full string serializer would require another document-sized allocation.
    # The supported artifact writer must operate through JSON's stream encoder.
    def whole_string(*args, **kwargs):
        pytest.fail('whole-document string allocation')

    monkeypatch.setattr(artifacts.json, 'dumps', whole_string)
    path = tmp_path/'artifact.json'
    artifacts.write_json_atomic(path, {'rows': [list(range(50)) for _ in range(100)]})
    assert json.loads(path.read_text())['rows'] == [list(range(50)) for _ in range(100)]


@pytest.mark.parametrize('failure', ['serialization', 'write', 'replace'])
def test_partial_stream_does_not_replace_previous_artifact(tmp_path, monkeypatch, failure):
    path = tmp_path/'artifact.json'
    path.write_bytes(b'previous verified artifact\n')
    previous = path.read_bytes()
    document = {'data': list(range(50)), 'later': [float('nan')]}
    expected = ValueError
    if failure == 'write':
        expected = OSError

        def interrupted(data, stream, **settings):
            stream.write('{"partial":')
            raise OSError('output full')

        monkeypatch.setattr(artifacts.json, 'dump', interrupted)
    elif failure == 'replace':
        expected = OSError
        document['later'] = []

        def failed_replace(source, destination):
            assert Path(destination) == path
            raise OSError('rename interrupted')

        monkeypatch.setattr(artifacts.os, 'replace', failed_replace)
    with pytest.raises(expected):
        artifacts.write_json_atomic(path, document)
    assert path.read_bytes() == previous
    assert sorted(p.name for p in tmp_path.iterdir()) == ['artifact.json']
