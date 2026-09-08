Deze repository bevat code t.b.v. analyses voor de Prognosetafels gepubliceerd door het Actuarieel Genootschap.

**Clone de repo**
- Zorg dat Git geïnstalleerd zijn
- Kopieer de https-link
- Clone de repo in de bestemde folder door in de commandprompt: git clone {link}

**Werkwijze**
- Werk nooit direct in de main. Maak een aparte branch aan voor elke deelopdracht. Deze branches kunnen vervolgens gereviewed worden in Pull Requests. Zo is overzichtelijk te zien wat de wijzigingen zijn geweest en kan gemakkelijk met elkaar overlegd worden.
- Aan het begin van de script kan je een runname opgeven. Vervolgens worden de resultaten in een gerelateerde map in "Runs/" geoutput.


**Hoe werkt Git?**
Git is een systeem waarmee je wijzigingen in code kunt bijhouden en waarmee meerdere mensen samen aan een project kunnen werken. Vaak wordt Git gebruikt samen met een online platform zoals GitHub.

De belangrijkste begrippen:

Commit: een opgeslagen momentopname van je wijzigingen. Je geeft er meestal een korte beschrijving bij, bijvoorbeeld "Loginpagina toegevoegd".
Branch: een aparte versie/lijn van de code waarin je kunt werken zonder de hoofdversie direct aan te passen. Bijvoorbeeld een branch voor een nieuwe functie.
Push: je lokale commits naar een online repository sturen, bijvoorbeeld naar GitHub.
Pull: de nieuwste wijzigingen van de online repository naar je eigen computer ophalen.
Merge: de wijzigingen van de ene branch samenvoegen met een andere branch, bijvoorbeeld een feature-branch met de main-branch.
Pull Request (PR): een verzoek om de wijzigingen van jouw branch samen te voegen met een andere branch. Teamleden kunnen de code eerst bekijken en feedback geven voordat de wijzigingen worden gemerged.

Simpel voorbeeld van process
Branch maken → code aanpassen → commit → push → Pull Request → controleren → merge

Daarna staat de nieuwe code in de hoofdbranch (main). Andere ontwikkelaars kunnen met pull de nieuwste versie ophalen.
