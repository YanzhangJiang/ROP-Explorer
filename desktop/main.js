const { app, BrowserWindow, dialog } = require('electron');
const { spawn } = require('node:child_process');
const fs = require('node:fs');
const path = require('node:path');

const BACKEND_PORT = Number(process.env.ROP_DESKTOP_PORT || 18088);
const COMPILED_BACKEND_STARTUP_TIMEOUT_MS = 90000;
const SOURCE_BACKEND_STARTUP_TIMEOUT_MS = 300000;

let backendProcess = null;
let isQuitting = false;

function backendExecutableName() {
  return process.platform === 'win32' ? 'rop-explorer-backend.exe' : 'rop-explorer-backend';
}

function juliaExecutableName() {
  return process.platform === 'win32' ? 'julia.exe' : 'julia';
}

function backendRootDir() {
  if (app.isPackaged) {
    const packagedCandidates = [
      path.join(process.resourcesPath, 'backend'),
      path.join(process.resourcesPath, 'ROPExplorerBackend'),
    ];

    for (const candidate of packagedCandidates) {
      if (fs.existsSync(candidate)) {
        return candidate;
      }
    }

    return packagedCandidates[0];
  }

  const compiledCandidate = path.join(__dirname, '..', 'dist', 'ROPExplorerBackend');
  const runtimeBundleCandidate = path.join(__dirname, '..', 'dist', 'ROPExplorerDesktopBackend');

  if (fs.existsSync(path.join(compiledCandidate, 'bin', backendExecutableName()))) {
    return compiledCandidate;
  }
  if (fs.existsSync(path.join(runtimeBundleCandidate, 'julia', 'bin', juliaExecutableName()))) {
    return runtimeBundleCandidate;
  }

  return compiledCandidate;
}

function compiledBackendExecutablePath(rootDir = backendRootDir()) {
  return path.join(rootDir, 'bin', backendExecutableName());
}

function compiledBackendPublicDir(rootDir = backendRootDir()) {
  return path.join(rootDir, 'share', 'rop-explorer', 'public');
}

function bundledJuliaExecutablePath(rootDir = backendRootDir()) {
  return path.join(rootDir, 'julia', 'bin', juliaExecutableName());
}

function bundledWebappDir(rootDir = backendRootDir()) {
  return path.join(rootDir, 'webapp');
}

function bundledServerPath(rootDir = backendRootDir()) {
  return path.join(bundledWebappDir(rootDir), 'server.jl');
}

function bundledDepotDir(rootDir = backendRootDir()) {
  return path.join(rootDir, 'depot');
}

function sourceBackendPublicDir(rootDir = backendRootDir()) {
  return path.join(bundledWebappDir(rootDir), 'public');
}

function userJuliaDepotDir() {
  return path.join(app.getPath('userData'), 'julia-depot');
}

function backendUrl() {
  return `http://127.0.0.1:${BACKEND_PORT}/`;
}

function sleep(ms) {
  return new Promise((resolve) => setTimeout(resolve, ms));
}

async function waitForBackend(timeoutMs) {
  const deadline = Date.now() + timeoutMs;
  let lastError = null;

  while (Date.now() < deadline) {
    try {
      const response = await fetch(backendUrl(), { method: 'GET' });
      if (response.ok) {
        return;
      }
      lastError = new Error(`Backend returned HTTP ${response.status}`);
    } catch (error) {
      lastError = error;
    }

    await sleep(1000);
  }

  throw new Error(`Backend did not become ready in ${Math.round(timeoutMs / 1000)}s.${lastError ? ` Last error: ${lastError.message}` : ''}`);
}

function stopBackend() {
  if (!backendProcess || backendProcess.killed) {
    return;
  }

  backendProcess.kill('SIGTERM');
}

function backendLaunchSpec() {
  const rootDir = backendRootDir();
  const compiledBackendPath = compiledBackendExecutablePath(rootDir);

  if (fs.existsSync(compiledBackendPath)) {
    return {
      command: compiledBackendPath,
      args: [],
      timeoutMs: COMPILED_BACKEND_STARTUP_TIMEOUT_MS,
      env: {
        ROP_PUBLIC_DIR: compiledBackendPublicDir(rootDir),
      },
    };
  }

  const juliaPath = bundledJuliaExecutablePath(rootDir);
  const serverPath = bundledServerPath(rootDir);
  if (!fs.existsSync(juliaPath) || !fs.existsSync(serverPath)) {
    throw new Error(
      `Missing backend resources. Checked ${compiledBackendPath} and ${juliaPath}.`,
    );
  }

  const writableDepotPath = userJuliaDepotDir();
  const logsDir = path.join(writableDepotPath, 'logs');
  fs.mkdirSync(logsDir, { recursive: true });

  return {
    command: juliaPath,
    args: [
      '--startup-file=no',
      `--project=${bundledWebappDir(rootDir)}`,
      serverPath,
    ],
    timeoutMs: SOURCE_BACKEND_STARTUP_TIMEOUT_MS,
    env: {
      ROP_PUBLIC_DIR: sourceBackendPublicDir(rootDir),
      JULIA_DEPOT_PATH: `${writableDepotPath}${path.delimiter}${bundledDepotDir(rootDir)}`,
      JULIA_PKG_OFFLINE: 'true',
      JULIA_PKG_PRECOMPILE_AUTO: '1',
      JULIA_PKG_SERVER: '',
      JULIA_HISTORY: path.join(logsDir, 'repl_history.jl'),
    },
  };
}

async function startBackend() {
  const launchSpec = backendLaunchSpec();

  backendProcess = spawn(launchSpec.command, launchSpec.args, {
    env: {
      ...process.env,
      ROP_PORT: String(BACKEND_PORT),
      ...launchSpec.env,
    },
    stdio: ['ignore', 'pipe', 'pipe'],
  });

  let stdoutBuffer = '';
  let stderrBuffer = '';

  backendProcess.stdout.on('data', (chunk) => {
    const text = chunk.toString();
    stdoutBuffer += text;
    process.stdout.write(`[backend] ${text}`);
  });

  backendProcess.stderr.on('data', (chunk) => {
    const text = chunk.toString();
    stderrBuffer += text;
    process.stderr.write(`[backend] ${text}`);
  });

  backendProcess.once('error', (error) => {
    dialog.showErrorBox('ROP Explorer backend error', error.stack || String(error));
    app.quit();
  });

  backendProcess.once('exit', (code, signal) => {
    backendProcess = null;

    if (isQuitting || signal === 'SIGTERM' || code === 0) {
      return;
    }

    const details = [stdoutBuffer.trim(), stderrBuffer.trim()].filter(Boolean).join('\n\n');
    dialog.showErrorBox(
      'ROP Explorer backend exited unexpectedly',
      details || `Backend exited with code ${code ?? 'unknown'}${signal ? ` (signal: ${signal})` : ''}.`,
    );
    app.quit();
  });

  await waitForBackend(launchSpec.timeoutMs);
}

async function createMainWindow() {
  const window = new BrowserWindow({
    width: 1440,
    height: 900,
    show: false,
    autoHideMenuBar: true,
    webPreferences: {
      contextIsolation: true,
      sandbox: true,
    },
  });

  window.once('ready-to-show', () => {
    window.show();
  });

  await window.loadURL(backendUrl());
}

app.whenReady().then(async () => {
  try {
    await startBackend();
    await createMainWindow();
  } catch (error) {
    dialog.showErrorBox('Failed to start ROP Explorer', error.stack || String(error));
    app.quit();
  }
});

app.on('activate', async () => {
  if (BrowserWindow.getAllWindows().length === 0) {
    await createMainWindow();
  }
});

app.on('before-quit', () => {
  isQuitting = true;
  stopBackend();
});

app.on('window-all-closed', () => {
  if (process.platform !== 'darwin') {
    app.quit();
  }
});
