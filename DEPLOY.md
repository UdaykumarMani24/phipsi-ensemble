# Putting the app on GitHub and Android

## 1. Upload to GitHub (no command line needed)

1. Sign in at https://github.com and click **New repository**.
2. Name it `phipsi-ensemble`, choose **Public**, tick nothing else, click **Create repository**.
3. On the new page click **uploading an existing file**.
4. Unzip `phipsi-ensemble.zip` on your computer, open the folder, select **everything inside it**
   (index.html, css, js, icons, vendor, tests, README.md ...) and drag it into the browser.
   Upload the *contents*, not the folder itself, so `index.html` sits at the top level.
5. Click **Commit changes**.

## 2. Turn on the free website (GitHub Pages)

1. In the repository go to **Settings → Pages**.
2. Under *Build and deployment*: Source = **Deploy from a branch**, Branch = **main**, folder = **/ (root)**. Click **Save**.
3. After 1–2 minutes the app is live at
   `https://<your-username>.github.io/phipsi-ensemble/`
4. Put this link in the README, your report and your paper.

## 3. Use it as an Android app (free, immediate)

1. Open the link above in **Chrome on the Android phone**.
2. Tap the **Install app** button in the header (or Chrome menu ⋮ → **Install app** / **Add to Home screen**).
3. The app gets its own icon, opens full-screen, and works offline. Downloading structures still needs internet.

## 4. Publish on the Google Play Store (optional)

1. Go to https://www.pwabuilder.com, paste your GitHub Pages link, click **Start**.
2. Choose **Android → Generate package**. You get a signed `.aab` file, a signing key, and an `assetlinks.json` file.
   **Keep the signing key safe**; you need it for every future update.
3. Create a Google Play developer account (one-time fee) at https://play.google.com/console and upload the `.aab`.
4. Without `assetlinks.json` the app shows a small browser address bar at the top. To remove it, the file must be served at
   `https://<your-username>.github.io/.well-known/assetlinks.json`, which means creating a second repository named
   exactly `<your-username>.github.io` and putting it in a `.well-known` folder there (also add an empty `.nojekyll` file).
5. Play Store review usually takes a few days. You will need a privacy policy link; the app collects no data, so a short page saying so is enough.

## 5. Get a DOI for citing (Zenodo)

1. Sign in at https://zenodo.org with your GitHub account.
2. In Zenodo **GitHub** settings, switch on `phipsi-ensemble`.
3. On GitHub create a **Release** (e.g. `v1.0.0`). Zenodo archives it and gives a DOI. Add the DOI to `CITATION.cff` and the README.

## Updating later

Edit or re-upload files on GitHub; the website updates automatically. Installed apps pick up the change the next time they open online.
If you change the app's own files, also change `phipsi-v1` to `phipsi-v2` in `sw.js` so phones refresh their offline copy.
