"""Download Austrian ALBINA avalanche bulletin data"""
import urllib.request, json, os, re
from datetime import datetime, timedelta

os.makedirs('data/cryosphere/austrian_avalanche', exist_ok=True)

# The avalanche.report endpoint returned HTML - let's parse it
# Try CAAML XML format which EAWS uses
print("=== Trying ALBINA CAAML endpoints ===")
test_urls = [
    'https://avalanche.report/albina/api/bulletins?date=2024-01-15&lang=en',
    'https://avalanche.report/albina/api/bulletins?date=2024-01-15&format=json&lang=en',
    'https://avalanche.report/albina/api/bulletins/2024-01-15?lang=en',
    'https://avalanche.report/albina/api/bulletins/latest?lang=en',
    'https://avalanche.report/albina_files/2024-01-15/2024-01-15_en.json',
    'https://static.avalanche.report/albina_files/2024-01-15/2024-01-15_en.json',
    'https://avalanche.report/albina_files/latest/en.json',
]

for url in test_urls:
    try:
        req = urllib.request.Request(url, headers={
            'User-Agent': 'Mozilla/5.0',
            'Accept': 'application/json, application/xml, text/html'
        })
        with urllib.request.urlopen(req, timeout=10) as resp:
            content_type = resp.headers.get('Content-Type', '')
            data = resp.read()
            print("  {} -> {} bytes [{}]".format(url.split('report')[-1], len(data), content_type[:40]))
            text = data.decode('utf-8', errors='replace')
            
            # Check if JSON
            if 'json' in content_type or text.strip().startswith('{') or text.strip().startswith('['):
                try:
                    j = json.loads(text)
                    if isinstance(j, list):
                        print("    JSON list: {} items".format(len(j)))
                        if j and isinstance(j[0], dict):
                            print("    Keys: {}".format(list(j[0].keys())))
                    elif isinstance(j, dict):
                        print("    JSON dict keys: {}".format(list(j.keys())[:10]))
                except:
                    pass
            elif 'xml' in content_type or text.strip().startswith('<?xml'):
                print("    XML content, first 300 chars:")
                print("    {}".format(text[:300]))
            else:
                # Check for danger level info in HTML
                danger_matches = re.findall(r'danger[_-]?(?:level|rating)["\s:=]+(\d)', text, re.I)
                if danger_matches:
                    print("    Danger levels found: {}".format(danger_matches[:10]))
                
                # Check for any JSON embedded in the page
                json_match = re.search(r'<script[^>]*>.*?(?:bulletins?|forecast|danger)\s*[:=]\s*(\[.*?\]|\{.*?\})', text, re.S | re.I)
                if json_match:
                    print("    Embedded JSON found: {}".format(json_match.group(1)[:200]))
                    
                print("    HTML, first 500 chars:")
                print("    {}".format(text[:500]))
    except Exception as e:
        print("  {} -> {}".format(url.split('report')[-1], str(e)[:80]))

# Try the simple bulletin page and scrape
print("\n=== Trying to scrape bulletin page ===")
url = 'https://avalanche.report/albina/api/bulletins?date=2024-01-15&lang=en'
try:
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    with urllib.request.urlopen(req, timeout=15) as resp:
        html = resp.read().decode('utf-8', errors='replace')
        
        # Look for danger level patterns
        # ALBINA uses dangerRating with elevation info
        patterns = [
            r'"dangerRating"\s*:\s*"?(\d)"?',
            r'danger-level-(\d)',
            r'class="[^"]*danger[^"]*(\d)',
            r'rating["\s:]+(\d)',
            r'level["\s:]+(\d)',
        ]
        for pat in patterns:
            matches = re.findall(pat, html)
            if matches:
                print("Pattern '{}' found: {}".format(pat[:30], matches[:10]))
except Exception as e:
    print("Scrape error: {}".format(e))

# Try EAWS API directly
print("\n=== EAWS Bulletins API ===")
eaws_urls = [
    'https://api.avalanche.report/v2/bulletin/2024-01-15',
    'https://api.avalanche.report/bulletin/2024-01-15',
]
for url in eaws_urls:
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0', 'Accept': 'application/json'})
        with urllib.request.urlopen(req, timeout=10) as resp:
            data = resp.read()
            print("  {} -> {} bytes".format(url, len(data)))
    except Exception as e:
        print("  {} -> {}".format(url, str(e)[:80]))
