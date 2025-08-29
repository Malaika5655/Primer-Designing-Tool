import tkinter as tk
from tkinter import ttk, messagebox
from tkinter.scrolledtext import ScrolledText

# ------------------------
# SAME FUNCTIONS AS BEFORE
# ------------------------
def reverse_complement(seq: str) -> str:
    comp = str.maketrans("ATGC", "TACG")
    return seq.translate(comp)[::-1]

def gc_content(seq: str) -> float:
    if not seq: return 0.0
    g, c = seq.count('G'), seq.count('C')
    return (g + c) * 100.0 / len(seq)

def tm_wallace_or_long(seq: str) -> float:
    n = len(seq)
    if n == 0: return 0.0
    gc, at = seq.count('G') + seq.count('C'), seq.count('A') + seq.count('T')
    return (4 * gc + 2 * at) if n <= 14 else (64.9 + 41 * (gc - 16.4) / n)

def valid_dna(seq: str) -> bool:
    return all(b in "ATGC" for b in seq)

def find_primers_in_sequence(template, min_len, max_len, gc_min, gc_max, tm_min, tm_max):
    n, hits = len(template), []
    for i in range(n):
        for L in range(min_len, max_len + 1):
            j = i + L
            if j > n: break
            p = template[i:j]
            gc, tm = gc_content(p), tm_wallace_or_long(p)
            if gc_min <= gc <= gc_max and tm_min <= tm <= tm_max:
                hits.append({"seq": p, "start": i, "end": j, "len": L,
                             "gc": round(gc, 2), "tm": round(tm, 2)})
    return hits

def choose_best(primer_list, tm_min, tm_max, target_len=20):
    if not primer_list: return None
    tm_mid = (tm_min + tm_max) / 2.0
    return min(primer_list, key=lambda p: abs(p["tm"] - tm_mid) + 0.05*abs(p["len"] - target_len))

# ------------------------
# GUI Logic
# ------------------------
def run_design():
    seq = entry_sequence.get("1.0", tk.END).strip().upper()
    try:
        min_len, max_len = int(spin_min_len.get()), int(spin_max_len.get())
        gc_min, gc_max = float(spin_gc_min.get()), float(spin_gc_max.get())
        tm_min, tm_max = float(spin_tm_min.get()), float(spin_tm_max.get())
    except ValueError:
        messagebox.showerror("Error", "Invalid numeric constraints.")
        return

    if not seq: return messagebox.showwarning("Input Required", "Enter a DNA sequence.")
    if not valid_dna(seq): return messagebox.showerror("Invalid", "Use only A/T/G/C.")
    if min_len > max_len: return messagebox.showerror("Error", "Min length > Max length.")

    f_hits = find_primers_in_sequence(seq, min_len, max_len, gc_min, gc_max, tm_min, tm_max)
    rc = reverse_complement(seq)
    r_hits = find_primers_in_sequence(rc, min_len, max_len, gc_min, gc_max, tm_min, tm_max)

    f_best, r_best = choose_best(f_hits, tm_min, tm_max), choose_best(r_hits, tm_min, tm_max)

    result_best.config(state="normal"); result_best.delete("1.0", tk.END)
    result_list.config(state="normal"); result_list.delete("1.0", tk.END)

    if f_best:
        result_best.insert(tk.END, f"🟩 Forward Primer: {f_best['seq']} "
                                   f"(Len {f_best['len']}, GC {f_best['gc']}%, Tm {f_best['tm']}°C)\n")
    if r_best:
        result_best.insert(tk.END, f"🟪 Reverse Primer: {r_best['seq']} "
                                   f"(Len {r_best['len']}, GC {r_best['gc']}%, Tm {r_best['tm']}°C)\n")

    def show_top(title, items):
        result_list.insert(tk.END, f"{title}\n")
        for i, p in enumerate(items[:5], 1):
            result_list.insert(tk.END, f"  {i}. {p['seq']} | Len {p['len']} | GC {p['gc']}% | Tm {p['tm']}°C\n")
        result_list.insert(tk.END, "\n")

    show_top("Forward Candidates:", f_hits)
    show_top("Reverse Candidates:", r_hits)

    result_best.config(state="disabled")
    result_list.config(state="disabled")

def clear_all():
    entry_sequence.delete("1.0", tk.END)
    result_best.config(state="normal"); result_best.delete("1.0", tk.END); result_best.config(state="disabled")
    result_list.config(state="normal"); result_list.delete("1.0", tk.END); result_list.config(state="disabled")

# ------------------------
# GUI Layout (Compact Fit)
# ------------------------
root = tk.Tk()
root.title("Primer Designing Tool")
root.geometry("900x650")
root.configure(bg="#0b1021")

# Title
ttk.Label(root, text="🧬 Primer Designing Tool", font=("Segoe UI", 16, "bold"), foreground="#22D3EE", background="#0b1021").pack(pady=8)

# Input DNA
ttk.Label(root, text="Enter DNA Sequence:", foreground="white", background="#0b1021").pack()
entry_sequence = ScrolledText(root, height=4, wrap="word", font=("Consolas", 11),
                              bg="#0f172a", fg="white", insertbackground="white")
entry_sequence.pack(fill="x", padx=12, pady=6)

# Constraints
frame_opts = tk.Frame(root, bg="#111827"); frame_opts.pack(fill="x", padx=12, pady=6)
for lbl, var, dflt in [("Min Len", "spin_min_len", 18), ("Max Len", "spin_max_len", 25),
                       ("GC Min", "spin_gc_min", 40.0), ("GC Max", "spin_gc_max", 60.0),
                       ("Tm Min", "spin_tm_min", 55.0), ("Tm Max", "spin_tm_max", 65.0)]:
    ttk.Label(frame_opts, text=lbl+":", foreground="#9CA3AF", background="#111827").pack(side="left", padx=4)
    sb = ttk.Spinbox(frame_opts, from_=0, to=100, width=6); sb.set(dflt); sb.pack(side="left", padx=4)
    globals()[var] = sb

# Buttons
btn_frame = tk.Frame(root, bg="#0b1021"); btn_frame.pack(pady=8)
ttk.Button(btn_frame, text="Find Primers", command=run_design).pack(side="left", padx=6)
ttk.Button(btn_frame, text="Clear", command=clear_all).pack(side="left", padx=6)

# Best Results
ttk.Label(root, text="Best Primers:", foreground="#A78BFA", background="#0b1021").pack()
result_best = ScrolledText(root, height=5, font=("Consolas", 11), bg="#0f172a", fg="#22D3EE")
result_best.pack(fill="x", padx=12, pady=6)
result_best.config(state="disabled")

# Candidate List
ttk.Label(root, text="Top Candidates:", foreground="#FF6B6B", background="#0b1021").pack()
result_list = ScrolledText(root, height=10, font=("Consolas", 11), bg="#0f172a", fg="white")
result_list.pack(fill="both", expand=True, padx=12, pady=6)
result_list.config(state="disabled")

root.mainloop()
