import sys, re, inspect
sys.path.insert(0, r"C:\Projects\Claude4\penta")
import shapes
for name in ("Z", "P"):
    ch = shapes.get_chart(name)
    src = inspect.getsource(ch["fn"])
    open(f"lambdify_src_{name}.txt", "w").write(src)
    print(name, "source length", len(src), "lines", src.count("\n"))
    # look for divisions, floats, rationals, sqrt, function calls
    body = src.split("\n", 1)[1]
    print("  '/' occurrences:", body.count("/"))
    print("  float literals:", re.findall(r"\d+\.\d*|\d+e[-+]?\d+", body)[:10])
    print("  function calls:", sorted(set(re.findall(r"([A-Za-z_][A-Za-z_0-9]*)\(", body))))
    print("  '**' exponents:", sorted(set(re.findall(r"\*\*\s*\(?(-?[0-9./]+)", body))))
    ints = [int(t) for t in re.findall(r"(?<![\w.])(\d+)(?![\w.])", body)]
    print("  max |int literal|:", max(ints), " bits", max(ints).bit_length())
    print("  monomial divisors:", {k: v for k, v in ch["mono"].items()})
    print("  first 400 chars:", body[:400].replace("\n", " | "))
