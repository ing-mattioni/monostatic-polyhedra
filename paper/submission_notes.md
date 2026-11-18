# Preparazione della pubblicazione C-only

Il paper corrente è *A homogeneous mono-unstable polyhedron with fourteen vertices*.
Versione 1.0, Claudio Mattioni, claudio.mattioni@gmail.com.

Gli allegati pronti sono in output/zenodo:

- mono-unstable-14.pdf
- mono-unstable-14-v1.0.zip

I valori completi del modulo sono in output/zenodo/ZENODO_FIELDS_IT.md e nella
copia locale outreach/zenodo_form.md. Il DOI riservato resta 10.5281/zenodo.24286668.
La data di pubblicazione va impostata al giorno della prima disponibilità pubblica.

Il vecchio ZIP e il manoscritto A–C non devono essere allegati al nuovo deposito.
Il precedente tag locale v1.0 identifica il vecchio snapshot, non questi file.
L'archivio corrente usa un elenco esplicito di file e il proprio manifest
CHECKSUMS_C.sha256; non viene prodotto con git archive del vecchio tag.
Il manoscritto e il supplemento precedenti sono conservati in
paper/history/abc-94a1f1e.

La bibliografia comprende il precedente pubblico a 16 vertici e l'archivio
Humeniuk. I risultati A/B non vengono rivendicati come nuovi.
Il materiale storico rimane nel repository privato e non viene cancellato.

Per una futura sottomissione arXiv, il sorgente principale è paper/latex/main.tex,
con equilibrium_ledger.tex e V14_plate.pdf. Tectonic 0.17.0 compila il PDF;
una compilazione effettiva con pdfLaTeX non è stata ancora verificata.
Questa limitazione non impedisce il deposito del PDF compilato su Zenodo.

La verifica matematica usa python reproduce.py, senza dipendenze esterne.
La revisione umana resta un passo successivo: la prova non è stata certificata
da un esperto esterno. Le vecchie bozze email su A/B sono superate; per una
richiesta di revisione usare il nuovo paper C-only e il DOI dopo la pubblicazione.

Nessuna email viene inviata e nessun deposito viene pubblicato da questi script.
