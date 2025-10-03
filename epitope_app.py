# epitope_app.py
import streamlit as st
import pandas as pd
from io import StringIO
from Bio import SeqIO  # from biopython
import base64

# -------------------------------
# Simple Rule-based Epitope Predictor
# -------------------------------
def predict_epitopes(seq, window=9, threshold=0.3):
    hydrophilic = set("KRDEQN")
    epitopes = []
    seq = seq.upper()
    for i in range(len(seq) - window + 1):
        peptide = seq[i:i + window]
        score = sum(aa in hydrophilic for aa in peptide) / window
        if score >= threshold:
            epitopes.append({
                "peptide": peptide,
                "start": i + 1,
                "end": i + window,
                "score": round(score, 2)
            })
    return epitopes

def highlight_sequence_html(seq, epitopes):
    """
    Return an HTML string where predicted epitopes are wrapped with a <span>
    so they appear highlighted. Overlapping epitopes are handled left-to-right.
    """
    markers = [False] * len(seq)
    for e in epitopes:
        for pos in range(e["start"] - 1, e["end"]):
            markers[pos] = True

    out = []
    i = 0
    while i < len(seq):
        if markers[i]:
            j = i
            while j < len(seq) and markers[j]:
                j += 1
            part = seq[i:j]
            out.append(f'<span style="background-color: #ffef9f; padding:2px; border-radius:3px;">{part}</span>')
            i = j
        else:
            out.append(seq[i])
            i += 1
    return "".join(out)

def fasta_to_seq(uploaded_file):
    content = uploaded_file.getvalue().decode("utf-8")
    records = list(SeqIO.parse(StringIO(content), "fasta"))
    if not records:
        return None, "No FASTA records found."
    # If multiple records, take first and warn
    if len(records) > 1:
        return str(records[0].seq), f"Multiple records found — using first record: {records[0].id}"
    return str(records[0].seq), None

def df_to_csv_download_link(df, filename="epitopes.csv"):
    csv = df.to_csv(index=False)
    b64 = base64.b64encode(csv.encode()).decode()
    href = f'<a href="data:file/csv;base64,{b64}" download="{filename}">⬇️ Download CSV</a>'
    return href

# -------------------------------
# Streamlit Web App
# -------------------------------
st.set_page_config(page_title="Epitope Prediction", layout="wide")
st.title("🧬 Epitope Prediction App")
st.markdown("""
A simple rule-based B-cell epitope predictor (sliding window, hydrophilic scoring).
You can paste a sequence or upload a FASTA file.  
""")

col1, col2 = st.columns([2,1])

with col1:
    seq_input = st.text_area("Paste protein sequence (single-letter codes) or leave blank to upload FASTA:", height=160)

    uploaded = st.file_uploader("Upload FASTA file (optional)", type=["fa","fasta","txt"])
    fasta_msg = None
    if uploaded:
        seq_from_fasta, fasta_msg = fasta_to_seq(uploaded)
        if seq_from_fasta is None:
            st.error(fasta_msg)
        else:
            if not seq_input.strip():
                seq_input = seq_from_fasta
            st.success("FASTA loaded." if fasta_msg is None else fasta_msg)

    window = st.number_input("Window length (peptide length)", min_value=5, max_value=25, value=9, step=1)
    threshold = st.slider("Hydrophilic fraction threshold", min_value=0.0, max_value=1.0, value=0.3, step=0.05)

    predict_btn = st.button("Predict Epitopes")

with col2:
    st.markdown("### How scoring works")
    st.markdown("- Sliding window of length `window` (default 9).  \n- Hydrophilic residues counted: `K R D E Q N`  \n- Score = fraction of hydrophilic residues in the window.  \n- Windows with score >= threshold are returned.")

if predict_btn:
    if not seq_input or not seq_input.strip():
        st.warning("Please enter a sequence or upload a FASTA file.")
    else:
        seq = "".join(seq_input.split()).upper()
        epitopes = predict_epitopes(seq, window=window, threshold=threshold)
        if not epitopes:
            st.error("No epitopes predicted with the current settings.")
        else:
            st.success(f"Predicted {len(epitopes)} epitope(s).")
            df = pd.DataFrame(epitopes)
            st.table(df)
            st.markdown(df_to_csv_download_link(df), unsafe_allow_html=True)

            # Highlight and display sequence
            html = highlight_sequence_html(seq, epitopes)
            st.markdown("### Sequence (highlighted predicted epitope regions):")
            st.markdown(f'<div style="font-family:monospace; white-space:pre-wrap; font-size:14px;">{html}</div>', unsafe_allow_html=True)

# Example and footer
st.info("Example sequence you can use: `MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYR`")
st.caption("This is a toy predictor for demonstration and teaching purposes only.")
