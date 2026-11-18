# Ricerca di un poliedro mono-instabile a 13 vertici

18 novembre 2026. In questa ricerca non è emerso un esempio con 13 vertici e
un solo equilibrio instabile rispetto al baricentro del volume uniforme.
I file di questa cartella appartengono agli esperimenti successivi al preprint.

## Una famiglia esclusa esattamente

Siano `p_0,...,p_13` i punti in `examples/V14_integer.json`. Per ciascun indice
`d`, sia `Q_d` l'inviluppo convesso dei tredici punti con indice diverso da
`d`. Per ogni trasformazione affine invertibile dello spazio, il solido
omogeneo immagine di `Q_d` ha almeno due vertici strettamente instabili.
Il risultato riguarda le quattordici forme così definite e tutte le loro
immagini affini. La questione dell'esistenza per tredici punti arbitrari
resta aperta in questo lavoro.

Per ricontrollare il risultato, dalla radice del repository, con la sola
libreria standard di Python:

```sh
python -I -S search_vertices/v13_followup/check_affine.py
```

Il comando rilegge le coordinate, ricostruisce i poligoni di faccia, le
adiacenze e il baricentro volumetrico. Controlla poi i certificati compressi
`affine_*.proof.json.gz` per 14 eliminazioni e 13 scelte del vertice esentato
per eliminazione. Le etichette dei 13 punti residui seguono il loro ordine
nel file originale. Per l'eliminazione 3, due alberi sono forniti in un
secondo sistema di coordinate invertibile (`affine_3_frame1.proof.json.gz`).
Il verificatore è distinto dal generatore e non importa il risolutore lineare.

### Perché gli alberi coprono tutte le deformazioni

1. Se la parte lineare della trasformazione è A, il baricentro si trasforma
   con il solido. Il segno del prodotto scalare per il vicino j del vertice i
   è dunque il segno di (p_j−p_i)ᵀ G (p_i−c), con G=AᵀA definita positiva e
   c baricentro della forma residua. La traslazione si cancella.
2. Un vertice è strettamente instabile se questi prodotti sono tutti negativi
   sui suoi spigoli. Se ci fosse al più un vertice strettamente instabile,
   si potrebbe esentare un indice r e scegliere per ogni altro vertice i un
   vicino j con prodotto nonnegativo. L'uguaglianza è ammessa in questa
   ipotesi, quindi anche i casi al confine sono inclusi.
3. Ogni biforcazione dell'albero sceglie un vertice vincolato e prova **tutti**
   i suoi vicini. Un ciclo di archi scelti è impossibile. Infatti lungo
   i→j la distanza quadratica dal baricentro, nella metrica G, cresce almeno
   di (p_j−p_i)ᵀ G (p_j−p_i), strettamente positiva.
4. All'altra tipologia di foglia si associa una somma a pesi interi positivi
   delle matrici simmetriche dei prodotti scelti, ciascuna normalizzata per
   un fattore positivo. Il suo opposto è semidefinito positivo e diverso da
   zero. Il checker lo stabilisce ispezionando tutti i minori principali e la
   traccia della matrice 3×3. Il prodotto di traccia della somma con G è
   quindi strettamente negativo, in contraddizione con i prodotti nonnegativi.

La scelta di quali alberi e pesi provare usa SciPy/HiGHS in virgola mobile.
Il controllo finale verifica i pesi con interi, ricostruisce ogni ramo e
richiede la copertura completa. Un certificato con un peso alterato in
memoria è stato respinto. Come controllo preliminare separato, `diagonal.py`
copre con poligoni razionali il caso delle sole dilatazioni sui tre assi e
ritrova una metrica valida per l'esempio a 14 vertici integro.

## Tentativi di costruzione

| Famiglia | Partenze | Esito |
|---|---:|---|
| Punti liberi in 3D, massimizzazione del peggior margine angolare | 144 | 93 partenze hanno conservato almeno un candidato con 13 punti estremi nelle tre fasi. Il miglior margine registrato è −0,007485. Il suo esempio razionale ha 13 vertici instabili. |
| Punti liberi in 3D, undici testimoni protetti e un dodicesimo mirato | 192 | 77 partenze hanno registrato un candidato con 13 vertici estremi e due instabili nel calcolo numerico. Il candidato conservato, seed 234, è stato verificato esattamente. |
| Due bracci scambiati da una rotazione di mezzo giro, con torsioni e ramificazioni | 48 | 1.320.704 valutazioni. Tutte le 48 configurazioni finali razionalizzate sono state controllate. Il miglior baricentro lungo l'asse ha quota +0,005968 in unità di altezza iniziale 1, mentre la costruzione cercava una quota negativa. |
| Trasformazioni proiettive delle 14 forme residue e variazione della metrica | 84 | 2.040.360 valutazioni. Tutti gli 84 candidati razionalizzati hanno 13 vertici e sono non degeneri. Il conteggio esatto minimo è U=3. |

La prima versione numerica della ricerca con torsione perdeva precisione su
corpi quasi piani e poteva restituire baricentri esterni. I risultati riportati
sono quelli della ripetizione con la coordinata trasversale normalizzata e con
ricontrollo razionale delle quote di tutti i candidati finali.

Il candidato `focused_144_192_best.json` ha 13 vertici, 33 spigoli e 22 facce.
Il rapporto `focused_144_192_exact.json` riconosce due facce stabili, due
spigoli di sella e due vertici instabili, 0 e 4, tutti non degeneri. Si può
verificare con:

```sh
python -I -S verify/exact_verify.py search_vertices/v13_followup/focused_144_192_best.json
```

Per ripetere le campagne con le dipendenze scientifiche del repository:

```sh
python search_vertices/v13_followup/spatial.py --start 0 --count 144 --workers 6
python search_vertices/v13_followup/spatial.py --start 144 --count 192 --workers 5 --focus
python search_vertices/v13_followup/twisted.py --count 48
python search_vertices/v13_followup/projective.py
python search_vertices/v13_followup/affine.py --budget 40000 --workers 7
python search_vertices/v13_followup/affine.py --indices 3 --roots 10 11 --frame 1 --budget 25000 --workers 1
```
