// Offline support: the app's own files are cached so it opens without a network.
// Downloads from RCSB / AlphaFold still need the internet.
const CACHE = 'phipsi-v1';
const FILES = [
  './', './index.html', './how-it-works.html', './css/style.css', './js/core.js', './js/app.js',
  './vendor/3Dmol-min.js', './manifest.webmanifest', './icons/icon-192.png', './icons/icon-512.png',
];

self.addEventListener('install', e => {
  e.waitUntil(caches.open(CACHE).then(c => c.addAll(FILES)).then(() => self.skipWaiting()));
});
self.addEventListener('activate', e => {
  e.waitUntil(caches.keys().then(keys => Promise.all(keys.filter(k => k !== CACHE).map(k => caches.delete(k))))
    .then(() => self.clients.claim()));
});
self.addEventListener('fetch', e => {
  const url = new URL(e.request.url);
  if (url.origin !== location.origin) return;            // never cache protein data from other sites
  e.respondWith(
    fetch(e.request).then(r => {                          // network first, so updates show up
      const copy = r.clone(); caches.open(CACHE).then(c => c.put(e.request, copy)); return r;
    }).catch(() => caches.match(e.request))
  );
});
