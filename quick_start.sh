#!/bin/bash
# Quick Start Script pour LODLink
# Usage: ./quick_start.sh

echo "╔══════════════════════════════════════════════════════════╗"
echo "║         LODLink — Démarrage Rapide                      ║"
echo "╚══════════════════════════════════════════════════════════╝"
echo ""

# Vérifier que les données existent
if [ ! -d "data" ]; then
    echo "❌ Erreur: Le dossier 'data/' n'existe pas"
    echo "   Veuillez créer un dossier 'data/' avec vos fichiers:"
    echo "   - pedfile.pro"
    echo "   - map"
    echo "   - freq"
    echo "   - genotyping"
    exit 1
fi

if [ ! -f "data/pedfile.pro" ]; then
    echo "❌ Erreur: data/pedfile.pro n'existe pas"
    exit 1
fi

echo "✅ Vérification des fichiers de données..."
echo "   📁 Pedigree: $(ls -lh data/pedfile.pro | awk '{print $5}')"
echo "   📁 Carte génétique: $(ls -lh data/map 2>/dev/null | awk '{print $5}' || echo 'non trouvé')"
echo "   📁 Fréquences: $(ls -lh data/freq 2>/dev/null | awk '{print $5}' || echo 'non trouvé')"
echo "   📁 Génotypage: $(ls -lh data/genotyping 2>/dev/null | awk '{print $5}' || echo 'non trouvé')"
echo ""

echo "🚀 Lancement de l'analyse genome-wide..."
echo "   Paramètres:"
echo "   - Modèle: dominant"
echo "   - Thinning: 0.5 cM"
echo "   - Extension: 3.0 Mb"
echo "   - HTML: activé"
echo ""

# Lancer l'analyse
uv run lodlink --html --extend-region 3.0 --thin 0.5

# Vérifier que l'analyse s'est bien passée
if [ $? -eq 0 ]; then
    echo ""
    echo "╔══════════════════════════════════════════════════════════╗"
    echo "║              ✅ ANALYSE TERMINÉE !                       ║"
    echo "╚══════════════════════════════════════════════════════════╝"
    echo ""
    echo "📊 Résultats disponibles dans: results/"
    echo ""
    echo "🌐 Pour voir le rapport HTML interactif:"
    if [[ "$OSTYPE" == "darwin"* ]]; then
        echo "   open results/linkage_results_interactive.html"
        echo ""
        echo "Voulez-vous ouvrir le rapport maintenant? (o/n)"
        read -r response
        if [[ "$response" =~ ^[Oo]$ ]]; then
            open results/linkage_results_interactive.html
        fi
    elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
        echo "   xdg-open results/linkage_results_interactive.html"
    else
        echo "   Ouvrez results/linkage_results_interactive.html dans votre navigateur"
    fi
    echo ""
    echo "📁 Autres fichiers:"
    echo "   - Pedigrees HaploPainter: results/haplopainter_*.png"
    echo "   - Vue genome-wide: results/genome_wide_lod.png"
    echo "   - Tableau résumé: results/lod_results_summary.tsv"
else
    echo ""
    echo "❌ L'analyse a échoué. Vérifiez les erreurs ci-dessus."
    exit 1
fi
