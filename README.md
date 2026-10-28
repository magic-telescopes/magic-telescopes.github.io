# MAGIC public webpage

This is the public webpage for the MAGIC Telescopes.
Due to the wqorld-wide disruptions related to vulnerabilities discovered within the PHP codebase
we are not currently able to fully open our entire webserver to the public.

This is a temporary and limited solution for external individuals who want
to access information about how to join and collaborate with MAGIC members
in scientific activities with the MAGIC telescopes during the next planned observation cycle.

## Maintainers

To test the website locally,

1. clone the repo
2. enter the resulting directory
3. ``python3 -m http.server 8000``
4. Open http://localhost:8000/

If you find missing or incomplete information, please open an issue
in this repository.

### Local development with VS Code

> [!TIP]
> This is the recommended way to work!

You can edit and preview the website locally using [Visual Studio Code](https://code.visualstudio.com/) and the [Live Server extension](https://marketplace.visualstudio.com/items?itemName=ritwickdey.LiveServer). Live Server serves the static website and automatically refreshes the browser when you save changes.

#### Prerequisites

- [Visual Studio Code](https://code.visualstudio.com/)
- The Live Server extension, recommended through `.vscode/extensions.json`

Open the repository folder in VS Code (**File → Open Folder…**). If prompted, install the recommended extensions.

#### Start the local website

1. Open `index.html` in VS Code.
2. Click **Go Live** in the status bar, or right-click the HTML file and select **Open with Live Server**.
3. Open the local website at [http://127.0.0.1:5500/](http://127.0.0.1:5500/), if it does not open automatically.

Live Server serves files from the repository root. This is required for JavaScript features such as loading the shared navigation menu with `fetch()`. Opening HTML files directly through a `file://` URL will not work correctly.

#### Editing and previewing changes

Edit the HTML, CSS, JavaScript, and image files directly in VS Code and save your changes. Live Server automatically refreshes the browser when relevant files change. CSS changes may be applied without a full page reload.

The Live Server configuration is stored in `.vscode/settings.json`. The default development port is `5500`.

To stop the server, click **Port: 5500** or the corresponding Live Server control in the VS Code status bar.

Alternatively, you can use the Python HTTP server task defined in `.vscode/tasks.json`, which serves the website on port `8000`. Do not run both servers on the same port if their configurations are changed.

#### Optional: Development container

A Dev Container configuration is provided in `.devcontainer/devcontainer.json` for contributors who prefer a consistent, containerized development environment.

To use it:

1. Install [Docker Desktop](https://www.docker.com/products/docker-desktop/), or another supported container runtime.
2. Install the [Dev Containers extension](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers).
3. Open the repository and run **Dev Containers: Reopen in Container** from the Command Palette.
4. Start Live Server using **Go Live**.
5. Open the forwarded local port in your browser.

The Dev Container should include the Live Server extension and forward its configured port to the host.
