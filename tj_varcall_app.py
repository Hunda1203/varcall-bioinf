import subprocess
import threading
import time
import tkinter as tk
from tkinter import filedialog, messagebox
from tkinter import filedialog
from tkinter import ttk
from tkinter.scrolledtext import ScrolledText


timer_running = False
start_time = 0


def choose_folder():
    folder = filedialog.askdirectory()
    if folder:
        folder_var.set(folder)


def write_output(text):
    output_box.insert(tk.END, text)
    output_box.see(tk.END)


def format_elapsed(seconds):
    hours = int(seconds // 3600)
    minutes = int((seconds % 3600) // 60)
    secs = int(seconds % 60)
    return f"{hours:02d}:{minutes:02d}:{secs:02d}"


def update_timer():
    if timer_running:
        elapsed = time.time() - start_time
        elapsed_var.set("Elapsed time: " + format_elapsed(elapsed))
        root.after(1000, update_timer)


def run_command(command, folder):
    write_output(f"\n▶ {command}\n")

    process = subprocess.Popen(
        command,
        shell=True,
        cwd=folder,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT
    )

    output = ""

    for line in process.stdout:
        output += line
        write_output(line)

    process.wait()

    if process.returncode == 0:
        write_output("✓ Done\n")
        return True, output
    else:
        write_output(f"✗ Error, return code: {process.returncode}\n")
        return False, output


def save_report(results, comparison_outputs, elapsed_time):

    read_type = read_type_var.get()

    path = filedialog.asksaveasfilename(
        title="Save Report",
        defaultextension=".txt",
        filetypes=[("Text files", "*.txt")],
        initialfile=f"varcall_report_{read_type}_{time.strftime('%Y_%m_%d_%H_%M')}.txt"
    )

    if not path:
        return None

    with open(path, "w", encoding="utf-8") as report:

        report.write("VARIANT CALLING REPORT\n")
        report.write("=" * 50 + "\n\n")

        report.write(f"Runtime: {elapsed_time}\n")
        report.write(f"Read type: {read_type}\n")
        report.write(f"Folder: {folder_var.get()}\n\n")

        report.write("Pipeline Steps\n")
        report.write("-" * 50 + "\n")

        for name, success in results:
            report.write(
                f"[{'PASS' if success else 'FAIL'}] {name}\n"
            )

        if comparison_outputs:

            report.write("\nComparison Results\n")
            report.write("-" * 50 + "\n\n")

            for title, output in comparison_outputs:
                report.write(title + "\n")
                report.write(output.strip() + "\n\n")

    return path

def save_report_popup(results, comparison_outputs, elapsed_time):

    path = save_report(
        results,
        comparison_outputs,
        elapsed_time
    )

    if path:
        messagebox.showinfo(
            "Report saved",
            f"Report saved as:\n\n{path}"
        )


def show_summary(results, comparison_outputs, elapsed_time):
    popup = tk.Toplevel(root)
    popup.title("Pipeline Results")
    popup.transient(root)
    popup.grab_set()
    popup.configure(bg="#ECECEC")

    popup_width = 1000
    popup_height = 700

    root.update_idletasks()

    root_x = root.winfo_x()
    root_y = root.winfo_y()
    root_width = root.winfo_width()
    root_height = root.winfo_height()

    popup_x = root_x + (root_width - popup_width) // 2
    popup_y = root_y + (root_height - popup_height) // 2

    popup.geometry(f"{popup_width}x{popup_height}+{popup_x}+{popup_y}")

    title = tk.Label(
        popup,
        text="Pipeline overview",
        font=("Segoe UI", 18),
        bg="#ECECEC",
        fg="#111111"
    )
    title.pack(anchor="w", padx=25, pady=(20, 10))

    text_box = ScrolledText(
        popup,
        font=("Consolas", 11),
        bg="white",
        fg="#111111",
        wrap="word",
        relief="flat",
        padx=15,
        pady=15
    )
    text_box.pack(fill="both", expand=True, padx=25, pady=(0, 20))

    text = f"Total runtime: {elapsed_time}\n\n"
    text += "Steps:\n"

    for name, success in results:
        text += f"{'✓' if success else '✗'} {name}\n"

    if comparison_outputs:
        text += "\nComparison results:\n"
        for title, output in comparison_outputs:
            text += f"\n--- {title} ---\n"
            text += output.strip() + "\n"

    text_box.insert("1.0", text)
    text_box.config(state="disabled")

    button_frame = tk.Frame(popup, bg="#ECECEC")
    button_frame.pack(pady=(0, 20))

    save_button = ttk.Button(
        button_frame,
        text="Save Report",
        style="Small.TButton",
        width=15,
        command=lambda: save_report_popup(results, comparison_outputs, elapsed_time)
    )
    save_button.pack(side="left", padx=10)

    close_button = ttk.Button(
        button_frame,
        text="Close",
        style="Small.TButton",
        width=15,
        command=popup.destroy
    )
    close_button.pack(side="left", padx=10)


def run_pipeline():
    global timer_running, start_time

    folder = folder_var.get()
    read_type = read_type_var.get().strip()

    if not folder or not read_type:
        write_output("Choose folder and read type first.\n")
        return

    ref = read_type + ".fasta"
    reads = read_type + "_simulated_reads.fasta"
    sam = read_type + ".sam"
    bam = read_type + "_filtered_sorted.bam"
    mpileup = read_type + "_mpileup.txt"

    steps = []

    if run_minimap_var.get():
        steps.append(("Running minimap2 alignment", f"minimap2 -x map-ont -a --eqx {ref} {reads} > {sam}", "normal"))

    if run_varcaller_var.get():
        steps.append(("Running custom variant caller", f"printf '{read_type}\\n' | python3 tj_varcall.py", "normal"))

    if run_bam_var.get():
        steps.append(("Creating sorted BAM", f"samtools view -b -F 2308 {sam} | samtools sort -o {bam}", "normal"))
        steps.append(("Indexing BAM", f"samtools index {bam}", "normal"))

    if run_mpileup_var.get():
        steps.append(("Creating mpileup.txt file", f"samtools mpileup -f {ref} {bam} > {mpileup}", "normal"))

    if run_comparison_var.get():
        if comparer_var.get() == "custom":
            steps.append(("Comparing results with 'mutated.csv' file", f"printf '{read_type}\\n' | python3 tj_eval.py", "comparison_key"))

        elif comparer_var.get() == "mpileup":
            steps.append(("Comparing results with 'mpileup.txt' file", f"printf '{read_type}\\n' | python3 tj_compare_mpileup.py", "comparison_mpileup"))

        elif comparer_var.get() == "both":
            steps.append(("Comparing results with 'mutated.csv' file", f"printf '{read_type}\\n' | python3 tj_eval.py", "comparison_key"))
            steps.append(("Comparing results with 'mpileup.txt' file", f"printf '{read_type}\\n' | python3 tj_compare_mpileup.py", "comparison_mpileup"))

    if not steps:
        write_output("No steps selected.\n")
        return

    run_button.config(state="disabled")
    progress_bar["value"] = 0
    status_var.set("Starting pipeline...")
    output_box.delete("1.0", tk.END)

    start_time = time.time()
    timer_running = True
    elapsed_var.set("Elapsed time: 00:00:00")
    root.after(1000, update_timer)

    results = []
    comparison_outputs = []
    total_steps = len(steps)

    for index, (label, command, step_type) in enumerate(steps, start=1):
        status_var.set(label)
        progress_bar["value"] = ((index - 1) / total_steps) * 100

        success, output = run_command(command, folder)
        results.append((label, success))

        if step_type == "comparison_key":
            comparison_outputs.append(("Comparison with 'mutated.csv' file", output))

        if step_type == "comparison_mpileup":
            comparison_outputs.append(("Comparison with 'mpileup.txt' file", output))

        if not success:
            timer_running = False
            final_time = format_elapsed(time.time() - start_time)
            status_var.set("Pipeline stopped because of an error.")
            run_button.config(state="normal")
            show_summary(results, comparison_outputs, final_time)
            return

        progress_bar["value"] = (index / total_steps) * 100

    timer_running = False
    final_time = format_elapsed(time.time() - start_time)
    elapsed_var.set("Elapsed time: " + final_time)

    status_var.set("Pipeline finished successfully.")
    write_output("\nPipeline finished successfully.\n")
    run_button.config(state="normal")
    show_summary(results, comparison_outputs, final_time)


def start_pipeline_thread():
    thread = threading.Thread(target=run_pipeline)
    thread.daemon = True
    thread.start()


def select_all_steps():
    run_minimap_var.set(True)
    run_varcaller_var.set(True)
    run_bam_var.set(True)
    run_mpileup_var.set(True)
    run_comparison_var.set(True)
    comparer_var.set("both")


def clear_all_steps():
    run_minimap_var.set(False)
    run_varcaller_var.set(False)
    run_bam_var.set(False)
    run_mpileup_var.set(False)
    run_comparison_var.set(False)


root = tk.Tk()
root.title("Variant Calling Pipeline")
root.geometry("1500x950+100+100")
root.minsize(1300, 850)
root.configure(bg="#ECECEC")

style = ttk.Style()
style.theme_use("clam")

style.configure("Apple.TButton", font=("Segoe UI", 16), padding=18, anchor="center")
style.configure("Small.TButton", font=("Segoe UI", 12), padding=10, anchor="center")
style.configure("TEntry", font=("Segoe UI", 13), padding=8)
style.configure("TCombobox", font=("Segoe UI", 13), padding=8)
style.configure("Horizontal.TProgressbar", thickness=18)

container = tk.Frame(root, bg="#ECECEC")
container.pack(fill="both", expand=True)

canvas = tk.Canvas(container, bg="#ECECEC", highlightthickness=0)
scrollbar = ttk.Scrollbar(container, orient="vertical", command=canvas.yview)

scrollable_frame = tk.Frame(canvas, bg="#ECECEC")
canvas_window = canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")


def update_scroll_region(event=None):
    canvas.configure(scrollregion=canvas.bbox("all"))


def resize_scrollable_frame(event):
    canvas.itemconfig(canvas_window, width=event.width)


scrollable_frame.bind("<Configure>", update_scroll_region)
canvas.bind("<Configure>", resize_scrollable_frame)

canvas.configure(yscrollcommand=scrollbar.set)

canvas.pack(side="left", fill="both", expand=True)
scrollbar.pack(side="right", fill="y")


def mouse_wheel_windows(event):
    canvas.yview_scroll(-1 * int(event.delta / 120), "units")


def mouse_wheel_linux_up(event):
    canvas.yview_scroll(-1, "units")


def mouse_wheel_linux_down(event):
    canvas.yview_scroll(1, "units")


canvas.bind_all("<MouseWheel>", mouse_wheel_windows)
canvas.bind_all("<Button-4>", mouse_wheel_linux_up)
canvas.bind_all("<Button-5>", mouse_wheel_linux_down)

main = tk.Frame(scrollable_frame, bg="#ECECEC")
main.pack(fill="both", expand=True, padx=40, pady=35)

title = tk.Label(
    main,
    text="Variant Calling Pipeline",
    font=("Segoe UI", 30, "bold"),
    bg="#ECECEC",
    fg="#111111"
)
title.pack(anchor="w")

subtitle = tk.Label(
    main,
    text="Run minimap2 alignment, custom variant caller, samtools mpileup and comparison scripts",
    font=("Segoe UI", 13),
    bg="#ECECEC",
    fg="#555555"
)
subtitle.pack(anchor="w", pady=(5, 30))

folder_var = tk.StringVar()
read_type_var = tk.StringVar(value="ecoli")
status_var = tk.StringVar(value="Ready.")
elapsed_var = tk.StringVar(value="Elapsed time: 00:00:00")

run_minimap_var = tk.BooleanVar(value=True)
run_varcaller_var = tk.BooleanVar(value=True)
run_bam_var = tk.BooleanVar(value=True)
run_mpileup_var = tk.BooleanVar(value=True)
run_comparison_var = tk.BooleanVar(value=True)
comparer_var = tk.StringVar(value="both")

card = tk.Frame(main, bg="white")
card.pack(fill="x", pady=(0, 25))

tk.Label(card, text="Select the folder with the data", bg="white", font=("Segoe UI", 13, "bold")).grid(
    row=0, column=0, sticky="w", padx=25, pady=(25, 8)
)

folder_entry = ttk.Entry(card, textvariable=folder_var)
folder_entry.grid(row=1, column=0, padx=25, pady=(0, 25), sticky="we")

browse_button = ttk.Button(card, text="Browse", style="Small.TButton", command=choose_folder)
browse_button.grid(row=1, column=1, padx=(0, 25), pady=(0, 25))

tk.Label(card, text="Choose read type", bg="white", font=("Segoe UI", 13, "bold")).grid(
    row=2, column=0, sticky="w", padx=25, pady=(0, 8)
)

read_dropdown = ttk.Combobox(
    card,
    textvariable=read_type_var,
    values=["ecoli", "lambda"],
    state="readonly",
    width=28
)
read_dropdown.grid(row=3, column=0, padx=25, pady=(0, 25), sticky="w")

card.columnconfigure(0, weight=1)

options_card = tk.Frame(main, bg="white")
options_card.pack(fill="x", pady=(0, 25))

options_card.grid_columnconfigure(0, weight=1)
options_card.grid_columnconfigure(1, weight=1)

tk.Label(
    options_card,
    text="Pipeline steps",
    bg="white",
    font=("Segoe UI", 14, "bold")
).grid(row=0, column=0, pady=(20, 10))

tk.Checkbutton(
    options_card,
    text="Run minimap2 alignment",
    variable=run_minimap_var,
    bg="white",
    font=("Segoe UI", 12)
).grid(row=1, column=0, sticky="w", padx=100)

tk.Checkbutton(
    options_card,
    text="Run variant caller",
    variable=run_varcaller_var,
    bg="white",
    font=("Segoe UI", 12)
).grid(row=2, column=0, sticky="w", padx=100)

tk.Checkbutton(
    options_card,
    text="Create sorted BAM + index",
    variable=run_bam_var,
    bg="white",
    font=("Segoe UI", 12)
).grid(row=3, column=0, sticky="w", padx=100)

tk.Checkbutton(
    options_card,
    text="Create 'mpileup.txt' file",
    variable=run_mpileup_var,
    bg="white",
    font=("Segoe UI", 12)
).grid(row=4, column=0, sticky="w", padx=100, pady=(0, 15))

ttk.Button(
    options_card,
    text="Select all",
    style="Small.TButton",
    command=select_all_steps
).grid(row=5, column=0, sticky="w", padx=80, pady=(0, 20))

ttk.Button(
    options_card,
    text="Clear all",
    style="Small.TButton",
    command=clear_all_steps
).grid(row=5, column=0, sticky="w", padx=220, pady=(0, 20))

tk.Label(
    options_card,
    text="Comparison method",
    bg="white",
    font=("Segoe UI", 14, "bold")
).grid(row=0, column=1, pady=(20, 10))

tk.Checkbutton(
    options_card,
    text="Run comparison",
    variable=run_comparison_var,
    bg="white",
    font=("Segoe UI", 12)
).grid(row=1, column=1, sticky="w", padx=100)

tk.Radiobutton(
    options_card,
    text="compare with 'mutated.csv' file",
    variable=comparer_var,
    value="custom",
    bg="white",
    font=("Segoe UI", 12)
).grid(row=2, column=1, sticky="w", padx=100)

tk.Radiobutton(
    options_card,
    text="compare with 'mpileup.txt' file",
    variable=comparer_var,
    value="mpileup",
    bg="white",
    font=("Segoe UI", 12)
).grid(row=3, column=1, sticky="w", padx=100)

tk.Radiobutton(
    options_card,
    text="Both methods",
    variable=comparer_var,
    value="both",
    bg="white",
    font=("Segoe UI", 12)
).grid(row=4, column=1, sticky="w", padx=100)

run_button = ttk.Button(
    main,
    text="Run Selected",
    style="Apple.TButton",
    command=start_pipeline_thread
)
run_button.pack(fill="x", pady=(0, 25))

progress_bar = ttk.Progressbar(
    main,
    orient="horizontal",
    mode="determinate",
    maximum=100,
    style="Horizontal.TProgressbar"
)
progress_bar.pack(fill="x", pady=(0, 10))

elapsed_label = tk.Label(
    main,
    textvariable=elapsed_var,
    font=("Segoe UI", 12),
    bg="#ECECEC",
    fg="#333333"
)
elapsed_label.pack(anchor="w", pady=(0, 8))

status_label = tk.Label(
    main,
    textvariable=status_var,
    font=("Segoe UI", 12),
    bg="#ECECEC",
    fg="#333333"
)
status_label.pack(anchor="w", pady=(0, 20))

output_frame = tk.Frame(main, bg="#ECECEC")
output_frame.pack(fill="x", pady=(0, 20))

output_box = ScrolledText(
    output_frame,
    font=("Consolas", 11),
    bg="#1E1E1E",
    fg="#FFFFFF",
    insertbackground="white",
    relief="flat",
    padx=15,
    pady=15,
    height=16
)
output_box.pack(fill="x")

root.mainloop()