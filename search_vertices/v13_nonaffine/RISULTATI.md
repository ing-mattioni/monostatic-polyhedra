# Ricerca non affine di un poliedro mono-instabile a 13 vertici

18 novembre 2026. Si cerca un solido convesso omogeneo con 13 vertici e un solo equilibrio instabile rispetto al baricentro del **volume**. Le cinque campagne hanno prodotto esempi esatti con due vertici instabili. Per una famiglia continua di modifiche dell'esempio a 14 vertici si può dimostrare che due instabili sono inevitabili.

## Dodici punti fissi e un punto libero

**Proposizione.** Si scelgano dodici punti qualsiasi dei quattordici elencati in [V14_integer.json](../../examples/V14_integer.json), mantenendo le coordinate indicate. Per ogni punto reale T dello spazio, se l'inviluppo convesso dei dodici punti e di T ha tredici vertici, almeno due dei suoi vertici sono strettamente instabili per il solido omogeneo.

In particolare, cancellare uno dei quattordici vertici e spostare uno solo dei tredici rimasti non può produrre l'esempio cercato. Il punto spostato può anche attraversare le facce o uscire dall'inviluppo originale. Il [precedente calcolo affine](../v13_followup/RISULTATI.md) riguarda un'altra famiglia. L'applicazione di una trasformazione affine arbitraria *dopo* la sostituzione e il problema generale a tredici punti restano aperti in questa ricerca.

### Il baricentro dopo l'aggiunta di T

Sia B l'inviluppo dei dodici punti fissi, c₀ il suo baricentro, D = 6 vol(B) e C = 4D c₀. Si triangola ciascuna faccia di B orientandola verso l'esterno. Per un triangolo con vertici a,b,d si pongono n = (b−a) × (d−a), o = n·a e S = a+b+d. La quantità w = max(0, n·T−o) è sei volte il volume del tetraedro aggiunto su quel triangolo quando è visibile da T. I tetraedri visibili ricoprono il volume aggiunto e hanno interni disgiunti. Il loro baricentro è (T+S)/4. Ne seguono le formule

    D_T = D + Σ w
    c(T) = [C + Σ w(T+S)] / (4 D_T).

Se tutti i pesi sono nulli, T appartiene a B e l'inviluppo ha al più dodici vertici.

### La verifica nel cubo centrale

Si copre il cubo chiuso ‖T‖∞ ≤ 45.056 (quarantacinquemilacinquantasei unità delle coordinate originali). In ciascun parallelepipedo si racchiudono i pesi e il baricentro in intervalli con estremi arrotondati verso l'esterno mediante aritmetica intera. Il codice interseca due inclusioni del baricentro, una ottenuta sommando i momenti rispetto a c₀ e l'altra osservando che la media delle somme S, pesata con i w, è una combinazione convessa.

Per ogni coppia di punti distinti i,j si considera il prodotto eᵢⱼ = (pⱼ−pᵢ)·(pᵢ−c(T)). Se tutti i prodotti uscenti da i sono strettamente negativi, il piano passante per pᵢ e normale a pᵢ−c(T) tocca il solido soltanto in pᵢ. È quindi un vertice strettamente instabile. Se invece eᵢⱼ ≥ 0, la distanza quadratica dal baricentro cresce da i a j di 2eᵢⱼ + ‖pⱼ−pᵢ‖² > 0.

Nel grafo del parallelepipedo si conserva un arco i→j ogni volta che l'intervallo calcolato per eᵢⱼ ammette valori non negativi. Per ogni T nel parallelepipedo il grafo contiene quindi tutti gli archi effettivamente possibili. In una componente fortemente connessa senza archi uscenti si prenda, per quel T, un punto alla massima distanza dal baricentro. Nessun suo prodotto verso un punto distinto può essere non negativo. Se esistono due di queste componenti, si ottengono due vertici strettamente instabili. Per riconoscerle basta controllare, mediante chiusura transitiva, che nessun nodo sia raggiungibile da tutti i nodi del grafo. In questo argomento si possono supporre distinti i tredici punti, come richiede l'ipotesi della proposizione.

Le celle indecise vengono bisecate, ricoprendo anche il piano di taglio. Sulle 91 scelte dei dodici punti il calcolo ha visitato **395.531 celle**, incluse quelle interne all'albero. Delle foglie, **197.806** certificano due componenti e **5** certificano T ∈ B. **Zero** foglie sono rimaste indecise. I [conteggi per ciascuna coppia rimossa](replacement_57_45000_cube45056.json) accompagnano il [programma che rigenera la copertura](replacement.py).

### Tutto lo spazio esterno al cubo

I quattordici vettori dati hanno norma euclidea al più 1.024. Trasliamo l'origine in c₀. I dodici punti fissi v sono ora contenuti nella palla di raggio R = 2.048. Indichiamo con y il vettore di T in queste coordinate e con L la sua norma. Fuori dal cubo vale L ≥ 45.056 − 1.024 > 20R.

Sia A l'area della proiezione ortogonale di B sul piano perpendicolare a y. Le fibre parallele a y hanno lunghezza al più 2R, quindi vol(B) ≤ 2RA. Nelle coordinate traslate gli offset o dei triangoli sono non negativi e la loro somma è D, per decomposizione in tetraedri con apice nell'origine. Ne segue D_T ≥ Σ max(0, n·y) = 2LA. L'ultima uguaglianza somma le aree proiettate dei triangoli rivolti verso y, che ricoprono una volta la proiezione a meno dei bordi. Dunque il volume ampliato è almeno LA/3.

La frazione di volume aggiunta η è almeno 1−6R/L. Per β = η/4 si ha allora 1/8 ≤ β ≤ 1/4. Poiché il baricentro di B è ora l'origine, il baricentro ampliato è c = β(y+s), dove s è una media pesata delle somme dei tre vertici dei triangoli visibili. In particolare ‖s‖ ≤ 3R. Per ogni punto fisso v,

    (y−v)·(v−c) ≤ −βL² + (1+4β)LR + 3βR²
                ≤ β(−L² + 12LR + 3R²) < 0.

Inoltre ‖y−c‖ ≥ (1−β)L−3βR, mentre ‖v−c‖ ≤ βL+(1+3β)R. La differenza è almeno L/2−5R/2 > 0. Il punto y è dunque il più lontano dal baricentro ed è strettamente instabile. Anche un punto fisso di distanza massima fra i dodici è strettamente instabile. Il prodotto verso ogni altro punto fisso è negativo per l'identità sulle distanze quadratiche, e quello verso y è negativo per la disuguaglianza appena dimostrata. Questo completa la copertura di tutto lo spazio.

## Ricerca con più vertici liberi

Il margine in tabella si calcola scegliendo un vertice r di distanza massima dal baricentro c. Per ogni i diverso da r si prende il massimo dei coseni tra pᵢ−c e le direzioni pⱼ−pᵢ, per tutti i j diversi da i. Il margine è il minimo di questi massimi. Un margine positivo, per un inviluppo con tredici veri vertici, garantisce che soltanto r sia instabile. A parità di due instabili, un margine negativo più vicino a zero indica che la migliore direzione disponibile al secondo instabile forma un angolo più vicino a 90° con il raggio uscente. Il valore zero richiede un esame delle degenerazioni. La separazione è la minima distanza fra due punti divisa per il diametro. I confronti dei coseni sono stati ordinati esattamente con il quadrato dotato di segno. I decimali sono arrotondati.

| Configurazione | Vertici instabili | Margine | Separazione |
|---|---:|---:|---:|
| Precedente candidato verificato con due instabili | 2 | −0,04927965 | 2,09506% |
| Nuovo candidato con separazione imposta ≥ 2,1% | 2 | −0,03563937 | 2,10031% |
| Nuovo candidato con separazione imposta ≥ 1,2% | 2 | −0,01968442 | 1,20164% |

Il miglior candidato della tabella ha 13 vertici, 33 spigoli, 22 facce, due equilibri stabili, due di sella e due instabili, tutti non degeneri. Gli instabili sono 0 e 12. Il vertice 12 è il più lontano dal baricentro e il vertice 1 offre la migliore direzione dal vertice 0. Rispetto al [punto di partenza](../v13_followup/spatial_0_144_best.json) sono cambiate esattamente le posizioni 0, 5, 7, 8, 9 e 10. Fra le sette posizioni fisse, 1, 2, 3 e 4 formano un tetraedro di volume non nullo. Ogni trasformazione affine che fissa questi quattro punti è l'identità, quindi una trasformazione affine comune non può spiegare lo spostamento dei sei vertici etichettati. Nel candidato a separazione maggiore i due instabili sono 0 e 1. Coordinate e controlli sono disponibili per [il miglior margine](moves_best.json), per [la separazione maggiore](moves_sep0.021_best.json) e nel [confronto esatto](audit.json).

| Esperimento | Avvii | Esito |
|---|---:|---|
| Rientro graduale di un vertice ausiliario e movimento indipendente dei tredici punti | 168 | Tre corpi iniziali a 14 vertici, tutte le cancellazioni e quattro modalità. I candidati sono valutati al baricentro dell'effettivo solido a 13 vertici. Miglior margine con due instabili −0,02843264. |
| Contorno a dodici angoli variabili e altezze determinate da programmazione lineare | 84 | 19.592 problemi lineari conteggiati. 62 avvii hanno dato coordinate finali a 13 vertici. Minimo di 12 vertici instabili fra questi esiti. |
| Riposizionamento successivo dei singoli vertici in tre dimensioni | 139 | 91 partenze dalla fusione di una coppia e 48 da forme precedenti. 138 esiti finali a 13 vertici. Miglior margine con due instabili −0,01968442. |
| Movimento coordinato di quattro vertici vicini, poi ottimizzazione di tutti e tredici con nuove scelte dei testimoni | 144 | 141 esiti finali a 13 vertici. Miglior margine con due instabili −0,02617460. |
| Ripetizione dei riposizionamenti con separazione minima del 2,1% per la selezione | 139 | 135 esiti finali a 13 vertici prima del filtro sulla separazione. Miglior margine ammesso con due instabili −0,03563937. |

I 674 avvii comprendono la ripetizione controllata degli stessi semi. Gli esiti di ottimizzazione sono evidenza numerica sulle forme esplorate. Il conteggio degli equilibri delle coordinate razionalizzate riportate in tabella è esatto. La possibilità di un esempio a tredici vertici di forma arbitraria rimane aperta.

## Riproduzione

Dalla radice del repository, con la sola libreria standard di Python, i comandi seguenti ricostruiscono rispettivamente la parte finita della prova, il confronto fra i candidati e la classificazione completa dell'esempio di miglior margine.

~~~sh
python -I -S search_vertices/v13_nonaffine/replacement.py --depth 57 --budget 45000 --cube 45056 --workers 8 --certify
python -I -S search_vertices/v13_nonaffine/audit.py
python -I -S verify/exact_verify.py search_vertices/v13_nonaffine/moves_best.json
~~~

Il primo comando termina con errore se resta anche una sola cella irrisolta. L'audit ha confrontato 364 volumi e baricentri con una ricostruzione razionale dell'inviluppo e ha controllato le esclusioni emesse sui parallelepipedi campionati. In altri 14 controlli, reinserendo ogni volta il vertice mancante nell'insieme degli altri tredici, il test non esclude erroneamente l'esempio mono-instabile noto.

I generatori numerici sono [continuation.py](continuation.py), [lift.py](lift.py), [move.py](move.py) e [reroute.py](reroute.py). I risultati per avvio sono raccolti nei file con suffisso `_runs.json` di questa cartella.
