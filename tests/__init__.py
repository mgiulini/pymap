"""Test directory."""
from pathlib import Path

tests_path = Path(__file__).resolve().parents[0]
reference_data = Path(tests_path, 'reference_data')
