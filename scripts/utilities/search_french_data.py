import urllib.request
import json

urls = [
    ('data.gouv.fr avalanche', 'https://www.data.gouv.fr/api/1/datasets/?q=avalanche&page_size=5'),
    ('data.gouv.fr EPA', 'https://www.data.gouv.fr/api/1/datasets/?q=EPA+avalanche&page_size=5'),
]

for name, url in urls:
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        resp = urllib.request.urlopen(req, timeout=15)
        data = json.loads(resp.read())
        if 'data' in data:
            print(name + ':')
            for item in data['data']:
                title = item.get('title', 'no title')
                org = item.get('organization', {})
                org_name = org.get('name', 'unknown') if org else 'unknown'
                print(f'  [{org_name}] {title}')
                resources = item.get('resources', [])
                for r in resources[:2]:
                    rtitle = r.get('title', '?')
                    rfmt = r.get('format', '?')
                    rurl = r.get('url', '?')[:120]
                    print(f'    -> {rtitle} ({rfmt}): {rurl}')
            print()
    except Exception as e:
        print(f'{name}: Error: {e}')
        print()
