#!/usr/bin/env python3
"""Fail-fast source/package audit for fastFGEE 0.3.0.9006."""
from __future__ import annotations

import json
import re
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def strip_r_strings_comments(text: str) -> str:
    out: list[str] = []
    i = 0
    quote: str | None = None
    escaped = False
    while i < len(text):
        c = text[i]
        if quote is not None:
            if escaped:
                escaped = False
                out.append(' ')
            elif c == '\\':
                escaped = True
                out.append(' ')
            elif c == quote:
                quote = None
                out.append(' ')
            else:
                out.append('\n' if c == '\n' else ' ')
            i += 1
            continue
        if c in ('"', "'", '`'):
            quote = c
            out.append(' ')
            i += 1
            continue
        if c == '#':
            while i < len(text) and text[i] != '\n':
                out.append(' ')
                i += 1
            continue
        out.append(c)
        i += 1
    if quote is not None:
        raise ValueError('unterminated R string/backtick')
    return ''.join(out)


def balanced(text: str, pairs: dict[str, str]) -> tuple[bool, str]:
    stack: list[tuple[str, int]] = []
    close_to_open = {v: k for k, v in pairs.items()}
    for i, c in enumerate(text):
        if c in pairs:
            stack.append((c, i))
        elif c in close_to_open:
            if not stack or stack[-1][0] != close_to_open[c]:
                return False, f'unmatched {c} at byte {i}'
            stack.pop()
    if stack:
        c, i = stack[-1]
        return False, f'unclosed {c} at byte {i}'
    return True, ''



def strip_cpp_strings_comments(text: str) -> str:
    """Remove C/C++ comments and string/character literals, preserving lines."""
    out: list[str] = []
    i = 0
    state = 'code'
    escaped = False
    while i < len(text):
        c = text[i]
        nxt = text[i + 1] if i + 1 < len(text) else ''
        if state == 'line_comment':
            if c == '\n':
                state = 'code'
                out.append('\n')
            else:
                out.append(' ')
            i += 1
            continue
        if state == 'block_comment':
            if c == '*' and nxt == '/':
                out.extend((' ', ' '))
                i += 2
                state = 'code'
            else:
                out.append('\n' if c == '\n' else ' ')
                i += 1
            continue
        if state in ('string', 'char'):
            if escaped:
                escaped = False
                out.append(' ')
            elif c == '\\':
                escaped = True
                out.append(' ')
            elif (state == 'string' and c == '"') or (state == 'char' and c == "'"):
                state = 'code'
                out.append(' ')
            else:
                out.append('\n' if c == '\n' else ' ')
            i += 1
            continue
        if c == '/' and nxt == '/':
            out.extend((' ', ' '))
            i += 2
            state = 'line_comment'
        elif c == '/' and nxt == '*':
            out.extend((' ', ' '))
            i += 2
            state = 'block_comment'
        elif c == '"':
            state = 'string'
            out.append(' ')
            i += 1
        elif c == "'":
            state = 'char'
            out.append(' ')
            i += 1
        else:
            out.append(c)
            i += 1
    if state in ('block_comment', 'string', 'char'):
        raise ValueError(f'unterminated C++ {state.replace("_", " ")}')
    return ''.join(out)


def cpp_rcpp_exports(path: Path) -> dict[str, int]:
    """Extract C++ function names/arities following Rcpp export attributes."""
    text = path.read_text()
    exports: dict[str, int] = {}
    marker_re = re.compile(r'//\s*\[\[Rcpp::export(?:\([^\]]*\))?\]\]')
    signature_re = re.compile(
        r'(?:[A-Za-z_][A-Za-z0-9_:<>]*[\s*&]+)+'
        r'([A-Za-z_][A-Za-z0-9_]*)\s*\((.*?)\)\s*\{',
        flags=re.S,
    )
    for marker in marker_re.finditer(text):
        match = signature_re.search(text, marker.end())
        if not match:
            raise ValueError(f'could not parse Rcpp export after byte {marker.start()}')
        name = match.group(1)
        args = split_top_level(match.group(2))
        if name in exports:
            raise ValueError(f'duplicate Rcpp export function: {name}')
        exports[name] = len(args)
    return exports

def rd_braces(text: str) -> tuple[bool, str]:
    depth = 0
    escaped = False
    for i, c in enumerate(text):
        if escaped:
            escaped = False
            continue
        if c == '\\':
            escaped = True
            continue
        if c == '{':
            depth += 1
        elif c == '}':
            depth -= 1
            if depth < 0:
                return False, f'unmatched }} at byte {i}'
    return (depth == 0, '' if depth == 0 else f'{depth} unclosed brace(s)')


def dcf_fields(path: Path) -> dict[str, str]:
    fields: dict[str, str] = {}
    key: str | None = None
    for line in path.read_text().splitlines():
        if line.startswith((' ', '\t')):
            if key is None:
                raise ValueError(f'continuation without field in {path}')
            fields[key] += ' ' + line.strip()
        elif ':' in line:
            key, value = line.split(':', 1)
            fields[key] = value.strip()
        elif line.strip():
            raise ValueError(f'malformed DCF line: {line}')
    return fields


def split_top_level(text: str, separator: str = ',') -> list[str]:
    """Split an R-like argument list without splitting nested expressions."""
    parts: list[str] = []
    start = 0
    stack: list[str] = []
    quote: str | None = None
    escaped = False
    pairs = {'(': ')', '[': ']', '{': '}'}
    for i, c in enumerate(text):
        if quote is not None:
            if escaped:
                escaped = False
            elif c == '\\':
                escaped = True
            elif c == quote:
                quote = None
            continue
        if c in ('"', "'", '`'):
            quote = c
            continue
        if c in pairs:
            stack.append(pairs[c])
            continue
        if stack and c == stack[-1]:
            stack.pop()
            continue
        if c == separator and not stack:
            parts.append(text[start:i].strip())
            start = i + 1
    parts.append(text[start:].strip())
    return [part for part in parts if part]


def formal_names(signature: str) -> list[str]:
    names: list[str] = []
    for part in split_top_level(signature):
        # The formal name is everything before the first top-level equals sign.
        lhs = split_top_level(part, '=')[0].strip()
        names.append(lhs)
    return names


def extract_macro(text: str, macro: str) -> str | None:
    marker = '\\' + macro + '{'
    start = text.find(marker)
    if start < 0:
        return None
    body_start = start + len(marker)
    depth = 1
    escaped = False
    for i in range(body_start, len(text)):
        c = text[i]
        if escaped:
            escaped = False
            continue
        if c == '\\':
            escaped = True
            continue
        if c == '{':
            depth += 1
        elif c == '}':
            depth -= 1
            if depth == 0:
                return text[body_start:i]
    return None


def extract_call_signature(text: str, function_name: str) -> str | None:
    match = re.search(rf'\b{re.escape(function_name)}\s*\(', text)
    if not match:
        return None
    start = match.end()
    depth = 1
    quote: str | None = None
    escaped = False
    for i in range(start, len(text)):
        c = text[i]
        if quote is not None:
            if escaped:
                escaped = False
            elif c == '\\':
                escaped = True
            elif c == quote:
                quote = None
            continue
        if c in ('"', "'", '`'):
            quote = c
        elif c == '(':
            depth += 1
        elif c == ')':
            depth -= 1
            if depth == 0:
                return text[start:i]
    return None


def main() -> None:
    checks: dict[str, object] = {}
    failures: list[str] = []

    desc = dcf_fields(ROOT / 'DESCRIPTION')
    checks['version'] = desc.get('Version')
    if desc.get('Version') != '0.3.0.9006':
        failures.append('DESCRIPTION version is not 0.3.0.9006')
    deps = ' '.join(desc.get(k, '') for k in ('Depends', 'Imports', 'Suggests', 'LinkingTo'))
    if re.search(r'\b(?:sanic|irregulAR1)\b', deps):
        failures.append('archived numerical dependency remains in DESCRIPTION fields')
    if 'Rcpp' not in desc.get('Imports', '') or desc.get('LinkingTo') != 'Rcpp':
        failures.append('Rcpp is not a direct Imports and LinkingTo dependency')
    if not re.search(r'refund\s*\(>=\s*0\.1-40\)', desc.get('Imports', '')):
        failures.append('refund minimum version is missing')

    # R lexical balance and top-level function definitions.
    definitions: dict[str, list[str]] = {}
    r_files = sorted((ROOT / 'R').glob('*.R'))
    syntax_r_files = r_files + sorted((ROOT / 'tests').rglob('*.R')) + \
        sorted((ROOT / 'inst/validation').glob('*.R'))
    for path in syntax_r_files:
        try:
            stripped = strip_r_strings_comments(path.read_text())
        except ValueError as e:
            failures.append(f'{path.relative_to(ROOT)}: {e}')
            continue
        ok, msg = balanced(stripped, {'(': ')', '[': ']', '{': '}'})
        if not ok:
            failures.append(f'{path.relative_to(ROOT)}: {msg}')
        if path.parent == ROOT / 'R':
            for m in re.finditer(r'^([A-Za-z.][A-Za-z0-9._]*)\s*<-\s*function\b',
                                 path.read_text(), flags=re.M):
                definitions.setdefault(m.group(1), []).append(str(path.relative_to(ROOT)))
    duplicate_defs = {k: v for k, v in definitions.items() if len(v) > 1}
    checks['top_level_function_definitions'] = len(definitions)
    checks['duplicate_top_level_definitions'] = duplicate_defs
    if duplicate_defs:
        failures.append(f'duplicate top-level function definitions: {sorted(duplicate_defs)}')
    for required in ('fgee', 'get_family_info', 'fgee_build_working_stats',
                     'fgee_update_working_cols_dt', '.fgee_fit_internal_core',
                     '.fgee_fit_internal'):
        if required not in definitions:
            failures.append(f'missing required function definition: {required}')
    if (ROOT / 'R/zzz_nuisance_integration.R').exists():
        failures.append('zzz_nuisance_integration.R still exists')

    # Public signature is read only up to its opening brace.
    public = (ROOT / 'R/fgee_public.R').read_text()
    m = re.search(r'^fgee\s*<-\s*function\s*\((.*?)\)\s*\{', public,
                  flags=re.M | re.S)
    if not m:
        failures.append('could not parse public fgee signature')
        signature = ''
        public_formals: list[str] = []
    else:
        signature = m.group(1)
        public_formals = formal_names(signature)
    checks['public_formals'] = public_formals
    prohibited = ['exact', 'gee.fit', 'max.iter', 'tune.method', 'working.engine']
    public_leaks = [x for x in prohibited if re.search(rf'(^|[\s,]){re.escape(x)}\s*=', signature)]
    checks['public_signature_prohibited_controls'] = public_leaks
    if public_leaks:
        failures.append(f'public fgee signature leaks controls: {public_leaks}')
    required_forcing = {
        'exact = FALSE', 'gee.fit = TRUE', 'max.iter = 1L',
        'tune.method = "one-step"', 'working.engine = "optimized"'
    }
    missing_forcing = sorted(x for x in required_forcing if x not in public)
    checks['missing_one_step_forcing'] = missing_forcing
    if missing_forcing:
        failures.append(f'public wrapper does not force one-step controls: {missing_forcing}')
    if 'identical(cv, "fastkfold")' not in public:
        failures.append('public wrapper can still trigger a legacy CV path')
    if 'out <- .fgee_fit_internal(' not in public:
        failures.append('public wrapper does not use the integrated nuisance path')
    if 'out <- .fgee_fit_internal_core(' in public:
        failures.append('public wrapper bypasses the integrated nuisance path')

    # Executable source cannot call the removed packages. Documentation may
    # mention them to explain removal and attribution.
    executable = '\n'.join(p.read_text() for p in r_files + sorted((ROOT / 'src').glob('*.[ch]pp')))
    removed_calls = re.findall(r'(?:sanic|irregulAR1)::[A-Za-z0-9._]+', executable)
    removed_guards = re.findall(r'requireNamespace\(["\'](?:sanic|irregulAR1)["\']', executable)
    checks['removed_dependency_calls'] = removed_calls
    checks['removed_dependency_guards'] = removed_guards
    if removed_calls or removed_guards:
        failures.append('executable source still references sanic or irregulAR1')
    stripped_executable = '\n'.join(
        strip_r_strings_comments(p.read_text()) for p in r_files
    )
    literal_guards = re.findall(r'if\s*\(\s*!?TRUE\s*\)', stripped_executable)
    checks['suspicious_literal_guards'] = literal_guards
    if literal_guards:
        failures.append(f'suspicious literal TRUE guards remain: {literal_guards}')
    if executable.count('fgee_iar1_apply_precision(') < 3:
        failures.append('internal irregular AR1 precision is not wired into all three legacy/optimized sites')

    # C++ lexical balance and Rcpp attribute-export arities.
    cpp_files = sorted((ROOT / 'src').glob('*.cpp'))
    cpp_exports: dict[str, int] = {}
    for path in cpp_files:
        try:
            stripped_cpp = strip_cpp_strings_comments(path.read_text())
            ok, msg = balanced(stripped_cpp, {'(': ')', '[': ']', '{': '}'})
            if not ok:
                failures.append(f'{path.relative_to(ROOT)}: {msg}')
            for name, arity in cpp_rcpp_exports(path).items():
                if name in cpp_exports:
                    failures.append(f'duplicate Rcpp export function: {name}')
                cpp_exports[name] = arity
        except ValueError as e:
            failures.append(f'{path.relative_to(ROOT)}: {e}')
    checks['cpp_attribute_exports'] = cpp_exports

    namespace = (ROOT / 'NAMESPACE').read_text()
    if 'sourceCpp' in namespace:
        failures.append('stale sourceCpp import remains')
    if 'importFrom(Rcpp,evalCpp)' not in namespace:
        failures.append('standard Rcpp evalCpp namespace import is missing')
    if 'useDynLib(fastFGEE, .registration = TRUE)' not in namespace:
        failures.append('registered DLL loading directive is missing')

    # Explicit Rcpp registration set.
    expected_arities = {
        '_fastFGEE_fgee_kron_inverse_kernel': 7,
        '_fastFGEE_fastk_fold_kernel': 9,
        '_fastFGEE_fgee_sympd_inverse_cpp': 1,
        '_fastFGEE_fgee_sympd_solve_cpp': 2,
        '_fastFGEE_fgee_iar1_precision_bands_cpp': 2,
        '_fastFGEE_fgee_iar1_precision_cpp': 2,
        '_fastFGEE_fgee_iar1_apply_precision_cpp': 3,
        '_fastFGEE_fgee_iar1_profile_nll_cpp': 3,
    }
    expected = list(expected_arities)
    rexp = (ROOT / 'R/RcppExports.R').read_text()
    cexp = (ROOT / 'src/RcppExports.cpp').read_text()
    missing_r = [x for x in expected if x not in rexp]
    missing_c = [x for x in expected if f'"{x}"' not in cexp]
    registered_pairs = re.findall(
        r'^\s*\{"(_fastFGEE_[^"]+)"\s*,[^,]+,\s*([0-9]+)\s*\}',
        cexp,
        flags=re.M,
    )
    registered_arities = {name: int(arity) for name, arity in registered_pairs}
    checks['expected_rcpp_entries'] = expected
    checks['expected_rcpp_arities'] = expected_arities
    checks['registered_rcpp_entries'] = list(registered_arities)
    checks['registered_rcpp_arities'] = registered_arities
    attribute_symbols = {
        f'_fastFGEE_{name}': arity for name, arity in cpp_exports.items()
    }
    checks['cpp_attribute_symbols'] = attribute_symbols
    if missing_r or missing_c or registered_arities != expected_arities:
        failures.append(
            'Rcpp registration mismatch; '
            f'missing R={missing_r}, missing C={missing_c}, '
            f'registered={registered_arities}'
        )
    if attribute_symbols != expected_arities:
        failures.append(
            'Rcpp attribute exports do not match the explicit expected set/arity: '
            f'{attribute_symbols}'
        )

    # Each R wrapper must expose the same number of arguments as its .Call.
    wrapper_arities: dict[str, int] = {}
    for symbol, expected_arity in expected_arities.items():
        function_name = symbol.removeprefix('_fastFGEE_')
        wrapper_match = re.search(
            rf'^(?:\.?{re.escape(function_name)})\s*<-\s*function\s*\((.*?)\)\s*\{{',
            rexp,
            flags=re.M | re.S,
        )
        if wrapper_match:
            wrapper_arities[symbol] = len(split_top_level(wrapper_match.group(1)))
        else:
            # The two historical wrappers intentionally have leading dots.
            wrapper_match = re.search(
                rf'^\.[A-Za-z0-9._]*{re.escape(function_name.split("_")[-1])}\s*<-\s*function\s*\((.*?)\)\s*\{{',
                rexp,
                flags=re.M | re.S,
            )
            if wrapper_match:
                wrapper_arities[symbol] = len(split_top_level(wrapper_match.group(1)))
    checks['R_wrapper_arities'] = wrapper_arities
    if any(wrapper_arities.get(sym) != arity for sym, arity in expected_arities.items()):
        failures.append(f'Rcpp R-wrapper arity mismatch: {wrapper_arities}')

    # Rd and vignette structural checks.
    for path in sorted((ROOT / 'man').glob('*.Rd')):
        ok, msg = rd_braces(path.read_text())
        if not ok:
            failures.append(f'{path.relative_to(ROOT)}: {msg}')
    fgee_rd = (ROOT / 'man/fgee.Rd').read_text()
    arguments_body = extract_macro(fgee_rd, 'arguments')
    usage_body = extract_macro(fgee_rd, 'usage')
    if arguments_body is None or extract_macro(fgee_rd, 'value') is None:
        failures.append('man/fgee.Rd is structurally incomplete')
    usage_signature = extract_call_signature(usage_body or '', 'fgee')
    if usage_signature is None:
        failures.append('man/fgee.Rd usage block could not be parsed')
        usage_formals: list[str] = []
    else:
        usage_formals = formal_names(usage_signature)
        leaked_rd = [x for x in prohibited if x in usage_formals]
        if leaked_rd:
            failures.append(f'man/fgee.Rd advertises removed controls: {leaked_rd}')
    argument_items = re.findall(r'\\item\{([^{}]+)\}', arguments_body or '')
    checks['rd_usage_formals'] = usage_formals
    checks['rd_argument_items'] = argument_items
    if public_formals and usage_formals != public_formals:
        failures.append(
            f'public/Rd usage formal mismatch: public={public_formals}, Rd={usage_formals}'
        )
    if public_formals and argument_items != public_formals:
        failures.append(
            f'public/Rd argument-item mismatch: public={public_formals}, Rd={argument_items}'
        )

    rmd = (ROOT / 'vignettes/fastFGEE.Rmd').read_text()
    fence_open = False
    fence_kind = None
    fence_count = 0
    r_chunk_count = 0
    for lineno, line in enumerate(rmd.splitlines(), start=1):
        stripped = line.strip()
        if not stripped.startswith('```'):
            continue
        fence_count += 1
        if not fence_open:
            fence_open = True
            fence_kind = stripped[3:].strip()
            if fence_kind.startswith('{r'):
                r_chunk_count += 1
        else:
            # Closing fences must be bare. A decorated fence while another is
            # open usually indicates a missing close earlier in the vignette.
            if stripped != '```':
                failures.append(
                    f'vignette opens a new decorated fence before closing the '
                    f'previous block at line {lineno}'
                )
            fence_open = False
            fence_kind = None
    checks['vignette_fence_count'] = fence_count
    checks['vignette_R_chunk_count'] = r_chunk_count
    if fence_open:
        failures.append(
            f'vignette has an unclosed fenced code block ({fence_kind!r})'
        )
    # Calls may discuss removed names in prose, but examples must not use them.
    code_chunks = '\n'.join(re.findall(r'```\{r[^}]*\}\n(.*?)```', rmd, flags=re.S))
    example_leaks = [x for x in prohibited if re.search(rf'\b{re.escape(x)}\s*=', code_chunks)]
    checks['vignette_example_prohibited_controls'] = example_leaks
    if example_leaks:
        failures.append(f'vignette examples use removed controls: {example_leaks}')

    # Artifacts and build ignore.
    stray = [str(p.relative_to(ROOT)) for p in ROOT.rglob('*')
             if p.is_file() and (p.name == 'Rplots.pdf' or p.suffix in ('.o', '.so', '.dll'))]
    checks['stray_build_artifacts'] = stray
    if stray:
        failures.append(f'stray build artifacts: {stray}')
    rb = (ROOT / '.Rbuildignore').read_text()
    if 'Rplots' not in rb:
        failures.append('.Rbuildignore does not exclude Rplots.pdf')

    checks['failures'] = failures
    checks['passed'] = not failures
    path = ROOT / 'validation/source_audit.json'
    path.write_text(json.dumps(checks, indent=2, sort_keys=True) + '\n')
    print(json.dumps(checks, indent=2, sort_keys=True))
    raise SystemExit(0 if not failures else 1)


if __name__ == '__main__':
    main()
