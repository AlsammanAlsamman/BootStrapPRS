"""bioconfigme

Centralized access to analysis and software configuration.

This module is intentionally small and stable: Snakemake rules should import
functions from here (via sys.path.append("utils")) and avoid opening YAML files
directly.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional, Sequence, Union

try:
	import yaml  # type: ignore
except Exception as exc:  # pragma: no cover
	raise ImportError(
		"PyYAML is required. Install it (e.g., `pip install PyYAML`)."
	) from exc


_REPO_ROOT = Path(__file__).resolve().parents[1]
_CONFIGS_DIR = _REPO_ROOT / "configs"
_ANALYSIS_PATH = _CONFIGS_DIR / "analysis.yml"
_SOFTWARE_PATH = _CONFIGS_DIR / "software.yml"


class _Missing:
	pass


_MISSING = _Missing()
_CACHE: Dict[str, Dict[str, Any]] = {}


def _read_text(path: Path) -> str:
	return path.read_text(encoding="utf-8")


def _clean_template_lines(text: str) -> str:
	"""Allow loading the verbatim templates even if they contain meta lines.

The required templates include human notes like "(Keep structure: ...)" which
are not valid YAML. We strip those lines during parsing while leaving the files
unchanged on disk.
"""

	cleaned: list[str] = []
	for line in text.splitlines():
		stripped = line.strip()
		if stripped.startswith("(") and stripped.endswith(")"):
			continue
		cleaned.append(line)
	return "\n".join(cleaned) + "\n"


def _load_yaml_mapping(path: Path) -> Dict[str, Any]:
	if not path.exists():
		raise FileNotFoundError(f"Missing config file: {path}")

	raw = _read_text(path)
	text = _clean_template_lines(raw)
	data = yaml.safe_load(text) or {}
	if not isinstance(data, dict):
		raise ValueError(f"Expected YAML mapping at top-level: {path}")
	return data


def load_analysis_config(*, reload: bool = False) -> Dict[str, Any]:
	"""Load configs/analysis.yml as a dict."""
	cache_key = "analysis"
	if reload or cache_key not in _CACHE:
		_CACHE[cache_key] = _load_yaml_mapping(_ANALYSIS_PATH)
	return _CACHE[cache_key]


def load_software_config(*, reload: bool = False) -> Dict[str, Any]:
	"""Load configs/software.yml as a dict."""
	cache_key = "software"
	if reload or cache_key not in _CACHE:
		_CACHE[cache_key] = _load_yaml_mapping(_SOFTWARE_PATH)
	return _CACHE[cache_key]


def get_results_dir() -> str:
	"""Return analysis.results_dir.

The value is returned exactly as specified in configs/analysis.yml.

If it is a relative path (e.g., "results"), Snakemake will interpret it
relative to the working directory where Snakemake is executed.
"""

	cfg = load_analysis_config()
	if "results_dir" not in cfg:
		raise KeyError("analysis.yml missing required key: results_dir")

	return str(cfg["results_dir"])


def get_software_module(tool: str) -> str:
	"""Return software[tool].module."""

	cfg = load_software_config()
	tool_cfg = cfg.get(tool)
	if not isinstance(tool_cfg, dict):
		raise KeyError(f"software.yml missing tool block: {tool}")
	module = tool_cfg.get("module")
	if module is None:
		raise KeyError(f"software.yml missing {tool}.module")
	return str(module)


def get_software_command(tool: str, default: Optional[str] = None) -> str:
	"""Return software[tool].command (or a sensible default).

	If the command field is missing, returns `default` if provided, otherwise
	returns `tool`.
	"""

	cfg = load_software_config()
	tool_cfg = cfg.get(tool)
	if not isinstance(tool_cfg, dict):
		raise KeyError(f"software.yml missing tool block: {tool}")
	command = tool_cfg.get("command")
	if command is None:
		return default if default is not None else tool
	return str(command)


def get_software_params(tool: str) -> Dict[str, Any]:
	"""Return software[tool].params as a dict (empty if missing)."""

	cfg = load_software_config()
	tool_cfg = cfg.get(tool)
	if not isinstance(tool_cfg, dict):
		raise KeyError(f"software.yml missing tool block: {tool}")
	params = tool_cfg.get("params")
	if params is None:
		return {}
	if not isinstance(params, dict):
		raise ValueError(f"software.yml {tool}.params must be a mapping")
	return params


def get_analysis_value(
	path: Union[str, Sequence[str]],
	default: Any = _MISSING,
) -> Any:
	"""Fetch a value from analysis.yml by path.

	Examples:
		get_analysis_value("harmonization.ref_fasta")
		get_analysis_value(["default_resources", "mem_mb"], default=32000)
"""

	if isinstance(path, str):
		keys = [k for k in path.split(".") if k]
	else:
		keys = list(path)

	cur: Any = load_analysis_config()
	for key in keys:
		if not isinstance(cur, dict) or key not in cur:
			if default is not _MISSING:
				return default
			raise KeyError(f"analysis.yml missing key path: {'.'.join(keys)}")
		cur = cur[key]
	return cur


