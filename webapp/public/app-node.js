// ROP Explorer — Node Edition Frontend (v3 — Dynamic Canvas)
const API = '';

// ===== State =====
const state = {
  sessionId: null,
  model: null,
  qK_syms: [],
};

// ===== Node Registry =====
let nodeIdCounter = 0;
const nodeRegistry = {}; // id → { type, el, data }
let connections = [];     // { fromNode, fromPort, toNode, toPort }

// ===== Port Types & Validation =====
const PORT_TYPES = {
  reactions: 'reactions',
  model: 'model',
  params: 'params',
  result: 'result',
};

const PORT_COLORS = {
  reactions: '#E8854A',
  model:    '#6C8CFF',
  params:   '#9B59B6',
  result:   '#17A2B8',
};

// ===== NODE_TYPES Registry =====
const NODE_TYPES = {
  'reaction-network': {
    category: 'input',
    headerClass: 'header-input',
    title: 'Reaction Network',
    inputs: [],
    outputs: [{ port: 'reactions', label: 'Reactions' }],
    defaultWidth: 280,
    createBody(nodeId) {
      return `
        <div class="reaction-header">
          <span class="reaction-header-label">Reaction</span>
          <span class="reaction-header-label reaction-header-kd">Kd (opt)</span>
          <span class="reaction-header-spacer"></span>
        </div>
        <div id="${nodeId}-reactions-list"></div>
        <button class="btn btn-small" onclick="addReactionRow('${nodeId}')">+ Add Reaction</button>
      `;
    },
    onInit(nodeId) {
      addReactionRow(nodeId, 'E + S <-> C_ES', 1e-3);
      addReactionRow(nodeId, 'E + P <-> C_EP', 1e-3);
    },
  },
  'model-builder': {
    category: 'process',
    headerClass: 'header-process',
    title: 'Model Builder',
    inputs: [{ port: 'reactions', label: 'Reactions' }],
    outputs: [{ port: 'model', label: 'Model' }],
    defaultWidth: 260,
    createBody(nodeId) {
      return `
        <div class="node-info" id="${nodeId}-model-info" style="display:none;">
          <pre id="${nodeId}-model-info-text"></pre>
        </div>
        <button class="btn btn-primary" onclick="buildModel('${nodeId}')">Build Model</button>
      `;
    },
    onInit(nodeId) {
      setupAutoModelBuild(nodeId);
    },
  },
  'model-summary': {
    category: 'result',
    headerClass: 'header-result',
    title: 'Model Summary',
    inputs: [{ port: 'model', label: 'Model' }],
    outputs: [],
    defaultWidth: 300,
    createBody(nodeId) {
      return `<div class="viewer-content" id="${nodeId}-content"><span class="text-dim">Connect to a Model Builder to see summary.</span></div>`;
    },
    async execute(nodeId) {
      const contentEl = document.getElementById(`${nodeId}-content`);
      if (!state.model) { contentEl.innerHTML = '<span class="text-dim">No model built yet.</span>'; return; }
      const m = state.model;
      contentEl.innerHTML = `
        <table>
          <tr><th>Property</th><th>Value</th></tr>
          <tr><td>Species (n)</td><td>${m.n}</td></tr>
          <tr><td>Totals (d)</td><td>${m.d}</td></tr>
          <tr><td>Reactions (r)</td><td>${m.r}</td></tr>
          <tr><td>Species</td><td>${m.x_sym.join(', ')}</td></tr>
          <tr><td>Totals</td><td>${m.q_sym.join(', ')}</td></tr>
          <tr><td>Constants</td><td>${m.K_sym.join(', ')}</td></tr>
        </table>
        <div style="margin-top:8px;"><strong>N matrix:</strong></div>
        <pre style="font-size:10px;color:#aaa;margin:4px 0;">${m.N.map(r => r.map(v => String(v).padStart(3)).join(' ')).join('\n')}</pre>
        <div><strong>L matrix:</strong></div>
        <pre style="font-size:10px;color:#aaa;margin:4px 0;">${m.L.map(r => r.map(v => String(v).padStart(3)).join(' ')).join('\n')}</pre>
      `;
    },
  },
  'vertices-table': {
    category: 'result',
    headerClass: 'header-result',
    title: 'Vertices Table',
    inputs: [{ port: 'model', label: 'Model' }],
    outputs: [],
    defaultWidth: 380,
    createBody(nodeId) {
      return `<div class="viewer-content" id="${nodeId}-content"><span class="text-dim">Waiting for model...</span></div>`;
    },
    async execute(nodeId) {
      const contentEl = document.getElementById(`${nodeId}-content`);
      setNodeLoading(nodeId, true);
      try {
        const data = await api('find_vertices', { session_id: state.sessionId });
        let html = '<table><thead><tr><th>#</th><th>Perm</th><th>Species</th><th>Type</th><th>Nullity</th></tr></thead><tbody>';
        data.vertices.forEach(v => {
          const typeTag = v.asymptotic
            ? '<span class="tag tag-asym">Asymp</span>'
            : '<span class="tag tag-nonasym">Non-A</span>';
          const singTag = v.singular
            ? ' <span class="tag tag-singular">Sing</span>'
            : ' <span class="tag tag-invertible">Inv</span>';
          const speciesStr = v.species ? v.species.join(', ') : '';
          html += `<tr><td>${v.idx}</td><td>[${v.perm.join(',')}]</td><td style="font-family:monospace;font-size:10px;">${speciesStr}</td><td>${typeTag}${singTag}</td><td>${v.nullity}</td></tr>`;
        });
        html += '</tbody></table>';
        contentEl.innerHTML = html;
      } catch (e) {
        contentEl.innerHTML = `<div class="node-error">${e.message}</div>`;
      }
      setNodeLoading(nodeId, false);
    },
  },
  'regime-graph': {
    category: 'result',
    headerClass: 'header-result',
    title: 'Regime Graph',
    inputs: [{ port: 'model', label: 'Model' }],
    outputs: [],
    defaultWidth: 420,
    createBody(nodeId) {
      return `<div class="viewer-content" id="${nodeId}-content"><span class="text-dim">Waiting for model...</span></div>`;
    },
    async execute(nodeId) {
      const contentEl = document.getElementById(`${nodeId}-content`);
      setNodeLoading(nodeId, true);
      try {
        const data = await api('build_graph', { session_id: state.sessionId });
        contentEl.innerHTML = `<div class="plot-container" id="${nodeId}-plot"></div>`;
        setTimeout(() => {
          plotRegimeGraph(data, `${nodeId}-plot`);
          setupPlotResize(nodeId, `${nodeId}-plot`);
        }, 50);
      } catch (e) {
        contentEl.innerHTML = `<div class="node-error">${e.message}</div>`;
      }
      setNodeLoading(nodeId, false);
    },
  },
  'siso-params': {
    category: 'parameter',
    headerClass: 'header-parameter',
    title: 'SISO Config',
    inputs: [{ port: 'model', label: 'Model' }],
    outputs: [{ port: 'params', label: 'Params' }],
    defaultWidth: 320,
    createBody(nodeId) {
      return `
        <div class="param-row">
          <label>Change qK:</label>
          <select id="${nodeId}-siso-select" class="auto-update"></select>
        </div>
        <div class="param-row">
          <label>Min (log10):</label>
          <input type="number" id="${nodeId}-min" value="-6" min="-20" max="20" step="0.5" class="auto-update">
        </div>
        <div class="param-row">
          <label>Max (log10):</label>
          <input type="number" id="${nodeId}-max" value="6" min="-20" max="20" step="0.5" class="auto-update">
        </div>
      `;
    },
    onInit(nodeId) {
      setupAutoUpdate(nodeId, 'siso-params');
    },
    async execute(nodeId) {
      // Populate the select with qK symbols
      const sel = document.getElementById(`${nodeId}-siso-select`);
      if (sel && state.qK_syms.length > 0) {
        const curVal = sel.value;
        sel.innerHTML = '';
        state.qK_syms.forEach(s => {
          const opt = document.createElement('option');
          opt.value = s; opt.textContent = s;
          sel.appendChild(opt);
        });
        if (curVal && state.qK_syms.includes(curVal)) sel.value = curVal;
        else sel.value = state.qK_syms[0];
      }

      // Store config in node data
      const info = nodeRegistry[nodeId];
      if (info) {
        info.data = info.data || {};
        info.data.config = {
          change_qK: sel ? sel.value : state.qK_syms[0],
          min: parseFloat(document.getElementById(`${nodeId}-min`)?.value || '-6'),
          max: parseFloat(document.getElementById(`${nodeId}-max`)?.value || '6')
        };
      }
    },
  },
  'siso-result': {
    category: 'result',
    headerClass: 'header-result',
    title: 'SISO Result',
    inputs: [{ port: 'params', label: 'Params' }],
    outputs: [],
    defaultWidth: 420,
    createBody(nodeId) {
      return `
        <button class="btn btn-small" onclick="computeSISOResult('${nodeId}')">Run</button>
        <div class="viewer-content" id="${nodeId}-content"><span class="text-dim">Click Run to compute SISO analysis</span></div>
      `;
    },
  },
  'siso-analysis': {
    category: 'viewer',
    headerClass: 'header-viewer',
    title: 'SISO Analysis',
    inputs: [{ port: 'model', label: 'Model' }],
    outputs: [],
    defaultWidth: 420,
    createBody(nodeId) {
      return `
        <div class="param-row">
          <label>Change qK:</label>
          <select id="${nodeId}-siso-select"></select>
        </div>
        <button class="btn btn-small" onclick="recomputeSISO('${nodeId}')">Recompute</button>
        <div class="viewer-content" id="${nodeId}-content"><span class="text-dim">Waiting for model...</span></div>
      `;
    },
    async execute(nodeId) {
      // Populate the select
      const sel = document.getElementById(`${nodeId}-siso-select`);
      if (sel && state.qK_syms.length > 0) {
        const curVal = sel.value;
        sel.innerHTML = '';
        state.qK_syms.forEach(s => {
          const opt = document.createElement('option');
          opt.value = s; opt.textContent = s;
          sel.appendChild(opt);
        });
        if (curVal && state.qK_syms.includes(curVal)) sel.value = curVal;
      }
      const changeQK = sel ? sel.value : state.qK_syms[0];
      if (!changeQK) return;

      const contentEl = document.getElementById(`${nodeId}-content`);
      setNodeLoading(nodeId, true);
      try {
        const data = await api('siso_paths', { session_id: state.sessionId, change_qK: changeQK });
        let html = `<div style="margin-bottom:8px;"><strong>${data.n_paths}</strong> paths, <strong>${data.sources.length}</strong> sources, <strong>${data.sinks.length}</strong> sinks</div>`;
        html += '<div class="path-list">';
        data.paths.forEach(p => {
          const permStr = p.perms.map(pr => `[${pr.join(',')}]`).join(' → ');
          html += `<div class="path-item" data-idx="${p.idx}" data-qk="${changeQK}" data-node="${nodeId}" onclick="selectSISOPath(this)">#${p.idx}: ${permStr}</div>`;
        });
        html += '</div>';
        html += `<div class="plot-container" id="${nodeId}-traj-plot" style="display:none;"></div>`;
        contentEl.innerHTML = html;
      } catch (e) {
        contentEl.innerHTML = `<div class="node-error">${e.message}</div>`;
      }
      setNodeLoading(nodeId, false);
    },
  },
  'rop-cloud': {
    category: 'viewer',
    headerClass: 'header-viewer',
    title: 'ROP Point Cloud',
    inputs: [{ port: 'reactions', label: 'Reactions' }, { port: 'model', label: 'Model' }],
    outputs: [],
    defaultWidth: 420,
    createBody(nodeId) {
      return `
        <div class="param-row">
          <label>Mode:</label>
          <select id="${nodeId}-sampling-mode" onchange="updateROPCloudMode('${nodeId}')">
            <option value="x_space">x-space closed-form</option>
            <option value="qk">qK sampling (legacy)</option>
          </select>
        </div>
        <div class="param-row">
          <label>Samples:</label>
          <input type="number" id="${nodeId}-samples" value="10000" min="100" max="100000" step="1000">
        </div>
        <div id="${nodeId}-xspace-params">
          <div class="param-row">
            <label>Target:</label>
            <select id="${nodeId}-target-species"></select>
          </div>
          <div class="param-row">
            <label>log10(x) min:</label>
            <input type="number" id="${nodeId}-logx-min" value="-6" min="-20" max="20" step="0.5">
          </div>
          <div class="param-row">
            <label>log10(x) max:</label>
            <input type="number" id="${nodeId}-logx-max" value="6" min="-20" max="20" step="0.5">
          </div>
        </div>
        <div id="${nodeId}-qk-params" style="display:none;">
          <div class="param-row">
            <label>Span:</label>
            <input type="number" id="${nodeId}-span" value="6" min="1" max="20">
          </div>
        </div>
        <button class="btn btn-small" onclick="recomputeROPCloud('${nodeId}')">Recompute</button>
        <div class="viewer-content" id="${nodeId}-content"><span class="text-dim">Waiting for input...</span></div>
      `;
    },
    onInit(nodeId) {
      updateROPCloudMode(nodeId);
    },
    async execute(nodeId) {
      const nSamples = parseInt(document.getElementById(`${nodeId}-samples`)?.value || '10000');
      const contentEl = document.getElementById(`${nodeId}-content`);
      const mode = document.getElementById(`${nodeId}-sampling-mode`)?.value || 'x_space';
      updateROPCloudMode(nodeId);
      setNodeLoading(nodeId, true);
      try {
        let data;
        if (mode === 'qk') {
          if (!state.sessionId) throw new Error('Build a model first, or switch to x-space mode');
          const modelConn = connections.find(c => c.toNode === nodeId && c.toPort === 'model');
          if (!modelConn) throw new Error('qK mode requires Model input connection');
          const span = parseInt(document.getElementById(`${nodeId}-span`)?.value || '6');
          data = await api('rop_cloud', {
            sampling_mode: 'qk',
            session_id: state.sessionId,
            n_samples: nSamples,
            span: span,
          });
        } else {
          const rxConn = connections.find(c => c.toNode === nodeId && c.toPort === 'reactions');
          if (!rxConn) throw new Error('x-space mode requires Reactions input connection');
          const { reactions } = getReactionsFromNode(rxConn.fromNode);
          if (!reactions.length) throw new Error('Add at least one reaction in the connected Reaction Network');
          refreshROPCloudTargetOptions(nodeId, reactions);
          const targetSpecies = document.getElementById(`${nodeId}-target-species`)?.value || '';
          const logxMin = parseFloat(document.getElementById(`${nodeId}-logx-min`)?.value || '-6');
          const logxMax = parseFloat(document.getElementById(`${nodeId}-logx-max`)?.value || '6');
          data = await api('rop_cloud', {
            sampling_mode: 'x_space',
            reactions: reactions,
            n_samples: nSamples,
            logx_min: logxMin,
            logx_max: logxMax,
            target_species: targetSpecies,
          });
        }
        contentEl.innerHTML = `<div class="plot-container" id="${nodeId}-plot"></div>`;
        setTimeout(() => {
          plotROPCloud(data, `${nodeId}-plot`);
          setupPlotResize(nodeId, `${nodeId}-plot`);
        }, 50);
      } catch (e) {
        contentEl.innerHTML = `<div class="node-error">${e.message}</div>`;
      }
      setNodeLoading(nodeId, false);
    },
  },
  'fret-heatmap': {
    category: 'viewer',
    headerClass: 'header-viewer',
    title: 'FRET Heatmap',
    inputs: [{ port: 'model', label: 'Model' }],
    outputs: [],
    defaultWidth: 420,
    createBody(nodeId) {
      return `
        <div class="param-row">
          <label>Grid size:</label>
          <input type="number" id="${nodeId}-grid" value="80" min="20" max="300">
        </div>
        <button class="btn btn-small" onclick="recomputeHeatmap('${nodeId}')">Recompute</button>
        <div class="viewer-content" id="${nodeId}-content"><span class="text-dim">Waiting for model (d=2 only)...</span></div>
      `;
    },
    async execute(nodeId) {
      const nGrid = parseInt(document.getElementById(`${nodeId}-grid`)?.value || '80');
      const contentEl = document.getElementById(`${nodeId}-content`);
      setNodeLoading(nodeId, true);
      try {
        const data = await api('fret_heatmap', {
          session_id: state.sessionId,
          n_grid: nGrid,
        });
        contentEl.innerHTML = `<div class="plot-container" id="${nodeId}-plot"></div>`;
        setTimeout(() => {
          plotHeatmap(data, `${nodeId}-plot`);
          setupPlotResize(nodeId, `${nodeId}-plot`);
        }, 50);
      } catch (e) {
        contentEl.innerHTML = `<div class="node-error">${e.message}</div>`;
      }
      setNodeLoading(nodeId, false);
    },
  },
  'parameter-scan-1d': {
    category: 'viewer',
    headerClass: 'header-viewer',
    title: 'Parameter Scan (1D)',
    inputs: [{ port: 'model', label: 'Model' }],
    outputs: [],
    defaultWidth: 420,
    createBody(nodeId) {
      return `
        <div class="param-row">
          <label>Scan parameter:</label>
          <select id="${nodeId}-param"></select>
        </div>
        <div class="param-row">
          <label>Range min:</label>
          <input type="number" id="${nodeId}-min" value="-6" step="0.5">
        </div>
        <div class="param-row">
          <label>Range max:</label>
          <input type="number" id="${nodeId}-max" value="6" step="0.5">
        </div>
        <div class="param-row">
          <label>Points:</label>
          <input type="number" id="${nodeId}-points" value="200" min="10" max="1000">
        </div>
        <div class="param-row">
          <label>Output expression:</label>
          <div style="display:flex;gap:4px;">
            <input type="text" id="${nodeId}-expr" placeholder="e.g., C_ES or 2*C_ES+E" style="flex:1;">
            <select id="${nodeId}-species-helper" onchange="insertSpecies1D('${nodeId}')" style="width:80px;">
              <option value="">Insert...</option>
            </select>
          </div>
        </div>
        <button class="btn btn-primary" onclick="runParameterScan1D('${nodeId}')">Run Scan</button>
        <div class="viewer-content" id="${nodeId}-content">
          <span class="text-dim">Connect to model and configure scan.</span>
        </div>
      `;
    },
    async execute(nodeId) {
      if (!state.model) return;

      const paramSelect = document.getElementById(`${nodeId}-param`);
      const speciesHelper = document.getElementById(`${nodeId}-species-helper`);

      if (paramSelect.options.length === 0) {
        state.model.q_sym.forEach(s => paramSelect.add(new Option(s, s)));
        state.model.K_sym.forEach(s => paramSelect.add(new Option(s, s)));
      }

      if (speciesHelper.options.length === 1) {
        state.model.x_sym.forEach(s => speciesHelper.add(new Option(s, s)));
      }
    },
  },
  'parameter-scan-2d': {
    category: 'viewer',
    headerClass: 'header-viewer',
    title: 'Parameter Scan (2D)',
    inputs: [{ port: 'model', label: 'Model' }],
    outputs: [],
    defaultWidth: 420,
    createBody(nodeId) {
      return `
        <div class="param-row">
          <label>X-axis parameter:</label>
          <select id="${nodeId}-param1"></select>
        </div>
        <div class="param-row">
          <label>X range:</label>
          <input type="number" id="${nodeId}-min1" value="-6" step="0.5" style="width:60px">
          to
          <input type="number" id="${nodeId}-max1" value="6" step="0.5" style="width:60px">
        </div>
        <div class="param-row">
          <label>Y-axis parameter:</label>
          <select id="${nodeId}-param2"></select>
        </div>
        <div class="param-row">
          <label>Y range:</label>
          <input type="number" id="${nodeId}-min2" value="-6" step="0.5" style="width:60px">
          to
          <input type="number" id="${nodeId}-max2" value="6" step="0.5" style="width:60px">
        </div>
        <div class="param-row">
          <label>Grid size:</label>
          <input type="number" id="${nodeId}-grid" value="80" min="20" max="200">
        </div>
        <div class="param-row">
          <label>Output expression:</label>
          <div style="display:flex;gap:4px;">
            <input type="text" id="${nodeId}-expr" placeholder="e.g., C_ES or 2*C_ES+E" style="flex:1;">
            <select id="${nodeId}-species-helper" onchange="insertSpecies2D('${nodeId}')" style="width:80px;">
              <option value="">Insert...</option>
            </select>
          </div>
        </div>
        <button class="btn btn-primary" onclick="runParameterScan2D('${nodeId}')">Run Scan</button>
        <div class="viewer-content" id="${nodeId}-content">
          <span class="text-dim">Connect to model and configure scan.</span>
        </div>
      `;
    },
    async execute(nodeId) {
      if (!state.model) return;

      const param1Select = document.getElementById(`${nodeId}-param1`);
      const param2Select = document.getElementById(`${nodeId}-param2`);
      const speciesHelper = document.getElementById(`${nodeId}-species-helper`);

      if (param1Select.options.length === 0) {
        [...state.model.q_sym, ...state.model.K_sym].forEach(s => {
          param1Select.add(new Option(s, s));
          param2Select.add(new Option(s, s));
        });
        if (state.model.q_sym.length >= 2) {
          param1Select.value = state.model.q_sym[0];
          param2Select.value = state.model.q_sym[1];
        }
      }

      if (speciesHelper.options.length === 1) {
        state.model.x_sym.forEach(s => speciesHelper.add(new Option(s, s)));
      }
    },
  },
  'rop-polyhedron': {
    category: 'viewer',
    headerClass: 'header-viewer',
    title: 'ROP Polyhedron',
    inputs: [{ port: 'model', label: 'Model' }],
    outputs: [],
    defaultWidth: 420,
    createBody(nodeId) {
      return `
        <div class="param-row">
          <label>Output expression:</label>
          <div style="display:flex;gap:4px;">
            <input type="text" id="${nodeId}-expr" placeholder="e.g., C_ES or 2*C_ES+E" style="flex:1;">
            <select id="${nodeId}-species-helper" onchange="insertSpeciesPoly('${nodeId}')" style="width:80px;">
              <option value="">Insert...</option>
            </select>
          </div>
        </div>
        <div class="param-row">
          <label>X-axis parameter:</label>
          <select id="${nodeId}-param1"></select>
        </div>
        <div class="param-row">
          <label>Y-axis parameter:</label>
          <select id="${nodeId}-param2"></select>
        </div>
        <div class="param-row">
          <label><input type="checkbox" id="${nodeId}-asymptotic" checked> Asymptotic only</label>
        </div>
        <div class="param-row">
          <label>Max vertices:</label>
          <input type="number" id="${nodeId}-max-vertices" value="1000" min="10" max="5000">
        </div>
        <button class="btn btn-primary" onclick="runROPPolyhedron('${nodeId}')">Compute Polyhedron</button>
        <div class="viewer-content" id="${nodeId}-content">
          <span class="text-dim">Connect to model and configure.</span>
        </div>
      `;
    },
    async execute(nodeId) {
      if (!state.model) return;

      const param1Select = document.getElementById(`${nodeId}-param1`);
      const param2Select = document.getElementById(`${nodeId}-param2`);
      const speciesHelper = document.getElementById(`${nodeId}-species-helper`);

      if (param1Select.options.length === 0) {
        [...state.model.q_sym, ...state.model.K_sym].forEach(s => {
          param1Select.add(new Option(s, s));
          param2Select.add(new Option(s, s));
        });
        if (state.model.q_sym.length >= 2) {
          param1Select.value = state.model.q_sym[0];
          param2Select.value = state.model.q_sym[1];
        }
      }

      if (speciesHelper.options.length === 1) {
        state.model.x_sym.forEach(s => speciesHelper.add(new Option(s, s)));
      }
    },
  },
  'scan-1d-params': {
    category: 'parameter',
    headerClass: 'header-parameter',
    title: 'Scan 1D Config',
    inputs: [{ port: 'model', label: 'Model' }],
    outputs: [{ port: 'params', label: 'Config' }],
    defaultWidth: 320,
    createBody(nodeId) {
      return `
        <div class="param-row">
          <label>Scan parameter:</label>
          <select id="${nodeId}-param" class="auto-update"></select>
        </div>
        <div class="param-row">
          <label>Range min:</label>
          <input type="number" id="${nodeId}-min" value="-6" step="0.5" class="auto-update">
        </div>
        <div class="param-row">
          <label>Range max:</label>
          <input type="number" id="${nodeId}-max" value="6" step="0.5" class="auto-update">
        </div>
        <div class="param-row">
          <label>Points:</label>
          <input type="number" id="${nodeId}-points" value="200" min="10" max="1000" class="auto-update">
        </div>
        <div class="param-row">
          <label>Output expression:</label>
          <div style="display:flex;gap:4px;">
            <input type="text" id="${nodeId}-expr" placeholder="e.g., C_ES or 2*C_ES+E" style="flex:1;" class="auto-update">
            <select id="${nodeId}-species-helper" onchange="insertSpecies1D('${nodeId}')" style="width:80px;">
              <option value="">Insert...</option>
            </select>
          </div>
        </div>
      `;
    },
    onInit(nodeId) {
      setupAutoUpdate(nodeId, 'scan-1d-params');
    },
    async execute(nodeId) {
      if (!state.model) return;

      const paramSelect = document.getElementById(`${nodeId}-param`);
      const speciesHelper = document.getElementById(`${nodeId}-species-helper`);

      if (paramSelect.options.length === 0) {
        state.model.q_sym.forEach(s => paramSelect.add(new Option(s, s)));
        state.model.K_sym.forEach(s => paramSelect.add(new Option(s, s)));
        // Save config after populating select
        triggerConfigUpdate(nodeId, 'scan-1d-params');
      }

      if (speciesHelper.options.length === 1) {
        state.model.x_sym.forEach(s => speciesHelper.add(new Option(s, s)));
      }
    },
  },
  'scan-1d-result': {
    category: 'result',
    headerClass: 'header-result',
    title: '1D Scan Result',
    inputs: [{ port: 'params', label: 'Config' }],
    outputs: [],
    defaultWidth: 420,
    createBody(nodeId) {
      return `
        <button class="btn btn-primary" onclick="executeScan1DResult('${nodeId}')">Run Scan</button>
        <div class="viewer-content" id="${nodeId}-content">
          <span class="text-dim">Connect to Scan 1D Config and click Run.</span>
        </div>
      `;
    },
    async execute(nodeId) {
      // This will be called by the button
    },
  },
  'rop-cloud-params': {
    category: 'parameter',
    headerClass: 'header-parameter',
    title: 'ROP Cloud Config',
    inputs: [{ port: 'reactions', label: 'Reactions' }, { port: 'model', label: 'Model' }],
    outputs: [{ port: 'params', label: 'Config' }],
    defaultWidth: 320,
    createBody(nodeId) {
      return `
        <div class="param-row">
          <label>Sampling mode:</label>
          <select id="${nodeId}-sampling-mode" onchange="updateROPCloudMode('${nodeId}')" class="auto-update">
            <option value="x_space">x-space closed-form</option>
            <option value="qk">qK sampling (legacy)</option>
          </select>
        </div>
        <div class="param-row">
          <label>Samples:</label>
          <input type="number" id="${nodeId}-samples" value="10000" min="100" max="100000" step="1000" class="auto-update">
        </div>
        <div id="${nodeId}-xspace-params">
          <div class="param-row">
            <label>Target:</label>
            <select id="${nodeId}-target-species" class="auto-update"></select>
          </div>
          <div class="param-row">
            <label>log10(x) min:</label>
            <input type="number" id="${nodeId}-logx-min" value="-6" min="-20" max="20" step="0.5" class="auto-update">
          </div>
          <div class="param-row">
            <label>log10(x) max:</label>
            <input type="number" id="${nodeId}-logx-max" value="6" min="-20" max="20" step="0.5" class="auto-update">
          </div>
        </div>
        <div id="${nodeId}-qk-params" style="display:none;">
          <div class="param-row">
            <label>Span:</label>
            <input type="number" id="${nodeId}-span" value="6" min="1" max="20" class="auto-update">
          </div>
        </div>
      `;
    },
    onInit(nodeId) {
      updateROPCloudMode(nodeId);
      setupAutoUpdate(nodeId, 'rop-cloud-params');
    },
    async execute(nodeId) {
      updateROPCloudMode(nodeId);
      const mode = document.getElementById(`${nodeId}-sampling-mode`)?.value || 'x_space';
      if (mode === 'x_space') {
        const rxConn = connections.find(c => c.toNode === nodeId && c.toPort === 'reactions');
        if (rxConn) {
          const { reactions } = getReactionsFromNode(rxConn.fromNode);
          if (reactions.length > 0) {
            refreshROPCloudTargetOptions(nodeId, reactions);
          }
        }
      }
    },
  },
  'rop-cloud-result': {
    category: 'result',
    headerClass: 'header-result',
    title: 'ROP Cloud Result',
    inputs: [{ port: 'params', label: 'Config' }],
    outputs: [],
    defaultWidth: 600,
    createBody(nodeId) {
      return `
        <button class="btn btn-primary" onclick="executeROPCloudResult('${nodeId}')">Run</button>
        <div class="viewer-content" id="${nodeId}-content">
          <span class="text-dim">Connect to ROP Cloud Config and click Run.</span>
        </div>
      `;
    },
    async execute(nodeId) {
      // This will be called by the button
    },
  },
  'fret-params': {
    category: 'parameter',
    headerClass: 'header-parameter',
    title: 'FRET Config',
    inputs: [{ port: 'model', label: 'Model' }],
    outputs: [{ port: 'params', label: 'Config' }],
    defaultWidth: 320,
    createBody(nodeId) {
      return `
        <div class="param-row">
          <label>Grid size:</label>
          <input type="number" id="${nodeId}-grid" value="80" min="20" max="300" class="auto-update">
        </div>
        <div class="param-row">
          <label>Min (log10):</label>
          <input type="number" id="${nodeId}-min" value="-6" min="-20" max="20" step="0.5" class="auto-update">
        </div>
        <div class="param-row">
          <label>Max (log10):</label>
          <input type="number" id="${nodeId}-max" value="6" min="-20" max="20" step="0.5" class="auto-update">
        </div>
      `;
    },
    onInit(nodeId) {
      setupAutoUpdate(nodeId, 'fret-params');
    },
    async execute(nodeId) {
      if (!state.model) return;
      // Store config in node data
      const info = nodeRegistry[nodeId];
      if (info) {
        info.data = info.data || {};
        info.data.config = {
          grid: parseInt(document.getElementById(`${nodeId}-grid`)?.value || '80'),
          min: parseFloat(document.getElementById(`${nodeId}-min`)?.value || '-6'),
          max: parseFloat(document.getElementById(`${nodeId}-max`)?.value || '6')
        };
      }
    },
  },
  'fret-result': {
    category: 'result',
    headerClass: 'header-result',
    title: 'FRET Result',
    inputs: [{ port: 'params', label: 'Config' }],
    outputs: [],
    defaultWidth: 600,
    createBody(nodeId) {
      return `
        <button class="btn btn-primary" onclick="executeFRETResult('${nodeId}')">Run</button>
        <div class="viewer-content" id="${nodeId}-content">
          <span class="text-dim">Connect to FRET Config and click Run.</span>
        </div>
      `;
    },
    async execute(nodeId) {
      // This will be called by the button
    },
  },
  'scan-2d-params': {
    category: 'parameter',
    headerClass: 'header-parameter',
    title: 'Scan 2D Config',
    inputs: [{ port: 'model', label: 'Model' }],
    outputs: [{ port: 'params', label: 'Config' }],
    defaultWidth: 320,
    createBody(nodeId) {
      return `
        <div class="param-row">
          <label>X-axis parameter:</label>
          <select id="${nodeId}-param1" class="auto-update"></select>
        </div>
        <div class="param-row">
          <label>X range min:</label>
          <input type="number" id="${nodeId}-min1" value="-6" step="0.5" class="auto-update">
        </div>
        <div class="param-row">
          <label>X range max:</label>
          <input type="number" id="${nodeId}-max1" value="6" step="0.5" class="auto-update">
        </div>
        <div class="param-row">
          <label>Y-axis parameter:</label>
          <select id="${nodeId}-param2" class="auto-update"></select>
        </div>
        <div class="param-row">
          <label>Y range min:</label>
          <input type="number" id="${nodeId}-min2" value="-6" step="0.5" class="auto-update">
        </div>
        <div class="param-row">
          <label>Y range max:</label>
          <input type="number" id="${nodeId}-max2" value="6" step="0.5" class="auto-update">
        </div>
        <div class="param-row">
          <label>Grid points:</label>
          <input type="number" id="${nodeId}-points" value="50" min="10" max="200" class="auto-update">
        </div>
        <div class="param-row">
          <label>Output expression:</label>
          <div style="display:flex;gap:4px;">
            <input type="text" id="${nodeId}-expr" placeholder="e.g., C_ES or 2*C_ES+E" style="flex:1;" class="auto-update">
            <select id="${nodeId}-species-helper" onchange="insertSpecies2D('${nodeId}')" style="width:80px;">
              <option value="">Insert...</option>
            </select>
          </div>
        </div>
      `;
    },
    onInit(nodeId) {
      setupAutoUpdate(nodeId, 'scan-2d-params');
    },
    async execute(nodeId) {
      if (!state.model) return;
      const param1Select = document.getElementById(`${nodeId}-param1`);
      const param2Select = document.getElementById(`${nodeId}-param2`);
      const speciesHelper = document.getElementById(`${nodeId}-species-helper`);
      if (param1Select.options.length === 0) {
        state.model.q_sym.forEach(s => {
          param1Select.add(new Option(s, s));
          param2Select.add(new Option(s, s));
        });
        state.model.K_sym.forEach(s => {
          param1Select.add(new Option(s, s));
          param2Select.add(new Option(s, s));
        });
        if (state.model.q_sym.length >= 2) {
          param1Select.value = state.model.q_sym[0];
          param2Select.value = state.model.q_sym[1];
        }
        // Save config after populating selects
        triggerConfigUpdate(nodeId, 'scan-2d-params');
      }
      if (speciesHelper.options.length === 1) {
        state.model.x_sym.forEach(s => speciesHelper.add(new Option(s, s)));
      }
    },
  },
  'scan-2d-result': {
    category: 'result',
    headerClass: 'header-result',
    title: '2D Scan Result',
    inputs: [{ port: 'params', label: 'Config' }],
    outputs: [],
    defaultWidth: 600,
    createBody(nodeId) {
      return `
        <button class="btn btn-primary" onclick="executeScan2DResult('${nodeId}')">Run Scan</button>
        <div class="viewer-content" id="${nodeId}-content">
          <span class="text-dim">Connect to Scan 2D Config and click Run.</span>
        </div>
      `;
    },
    async execute(nodeId) {
      // This will be called by the button
    },
  },
  'rop-poly-params': {
    category: 'parameter',
    headerClass: 'header-parameter',
    title: 'ROP Polyhedron Config',
    inputs: [{ port: 'model', label: 'Model' }],
    outputs: [{ port: 'params', label: 'Config' }],
    defaultWidth: 320,
    createBody(nodeId) {
      return `
        <div class="param-row">
          <label>Output expression:</label>
          <div style="display:flex;gap:4px;">
            <input type="text" id="${nodeId}-expr" placeholder="e.g., C_ES or 2*C_ES+E" style="flex:1;" class="auto-update">
            <select id="${nodeId}-species-helper" onchange="insertSpeciesPoly('${nodeId}')" style="width:80px;">
              <option value="">Insert...</option>
            </select>
          </div>
        </div>
        <div class="param-row">
          <label>X-axis parameter:</label>
          <select id="${nodeId}-param1" class="auto-update"></select>
        </div>
        <div class="param-row">
          <label>Y-axis parameter:</label>
          <select id="${nodeId}-param2" class="auto-update"></select>
        </div>
        <div class="param-row">
          <label><input type="checkbox" id="${nodeId}-asymptotic" checked class="auto-update"> Asymptotic only</label>
        </div>
        <div class="param-row">
          <label>Max vertices:</label>
          <input type="number" id="${nodeId}-max-vertices" value="1000" min="10" max="5000" class="auto-update">
        </div>
      `;
    },
    onInit(nodeId) {
      setupAutoUpdate(nodeId, 'rop-poly-params');
    },
    async execute(nodeId) {
      if (!state.model) return;
      const param1Select = document.getElementById(`${nodeId}-param1`);
      const param2Select = document.getElementById(`${nodeId}-param2`);
      const speciesHelper = document.getElementById(`${nodeId}-species-helper`);
      if (param1Select.options.length === 0) {
        [...state.model.q_sym, ...state.model.K_sym].forEach(s => {
          param1Select.add(new Option(s, s));
          param2Select.add(new Option(s, s));
        });
        if (state.model.q_sym.length >= 2) {
          param1Select.value = state.model.q_sym[0];
          param2Select.value = state.model.q_sym[1];
        }
      }
      if (speciesHelper.options.length === 1) {
        state.model.x_sym.forEach(s => speciesHelper.add(new Option(s, s)));
      }
    },
  },
  'rop-poly-result': {
    category: 'result',
    headerClass: 'header-result',
    title: 'ROP Polyhedron Result',
    inputs: [{ port: 'params', label: 'Config' }],
    outputs: [],
    defaultWidth: 600,
    createBody(nodeId) {
      return `
        <button class="btn btn-primary" onclick="executeROPPolyResult('${nodeId}')">Compute</button>
        <div class="viewer-content" id="${nodeId}-content">
          <span class="text-dim">Connect to ROP Polyhedron Config and click Compute.</span>
        </div>
      `;
    },
    async execute(nodeId) {
      // This will be called by the button
    },
  },
};

// Required predecessor chain for each node type
const PREREQ_CHAIN = {
  'model-builder': ['reaction-network'],
  'model-summary': ['reaction-network', 'model-builder'],
  'vertices-table': ['reaction-network', 'model-builder'],
  'regime-graph': ['reaction-network', 'model-builder'],
  'siso-analysis': ['reaction-network', 'model-builder'],
  'siso-params': ['reaction-network', 'model-builder'],
  'siso-result': ['siso-params'],
  'rop-cloud': ['reaction-network'],
  'fret-heatmap': ['reaction-network', 'model-builder'],
  'parameter-scan-1d': ['reaction-network', 'model-builder'],
  'parameter-scan-2d': ['reaction-network', 'model-builder'],
  'rop-polyhedron': ['reaction-network', 'model-builder'],
  'scan-1d-params': ['reaction-network', 'model-builder'],
  'scan-1d-result': ['scan-1d-params'],
  'rop-cloud-params': ['reaction-network'],
  'rop-cloud-result': ['rop-cloud-params'],
  'fret-params': ['reaction-network', 'model-builder'],
  'fret-result': ['fret-params'],
  'scan-2d-params': ['reaction-network', 'model-builder'],
  'scan-2d-result': ['scan-2d-params'],
  'rop-poly-params': ['reaction-network', 'model-builder'],
  'rop-poly-result': ['rop-poly-params'],
};

// ===== Canvas Interaction =====
const editor = document.getElementById('editor');
const canvas = document.getElementById('canvas');
const svgLayer = document.getElementById('svg-layer');

let panX = 0, panY = 0, isPanning = false, startPanX = 0, startPanY = 0;
let scale = 1.0;
const MIN_SCALE = 0.1;
const MAX_SCALE = 3.0;
let isDraggingNode = false, draggedNode = null, nodeOffsetX = 0, nodeOffsetY = 0;
let isWiring = false, wireStartSocket = null, wireStartIsOutput = true, tempWire = null;
let isResizing = false, resizeNode = null, resizeStartX = 0, resizeStartY = 0, resizeStartW = 0, resizeStartH = 0;

function applyViewportTransform() {
  // Keep a single transform source (canvas). svgLayer stays in canvas-local coordinates.
  canvas.style.transform = `translate(${panX}px, ${panY}px) scale(${scale})`;
  svgLayer.style.transform = 'none';
  editor.style.backgroundPosition = `${panX}px ${panY}px`;
  editor.style.backgroundSize = `${50 * scale}px ${50 * scale}px`;
}

// Canvas panning (middle/right mouse button, or left-click on blank area)
editor.addEventListener('mousedown', (e) => {
  if (e.button === 1 || e.button === 2) {
    isPanning = true;
    startPanX = e.clientX - panX;
    startPanY = e.clientY - panY;
    e.preventDefault();
  } else if (e.button === 0 && (e.target === editor || e.target === canvas || e.target === svgLayer)) {
    isPanning = true;
    startPanX = e.clientX - panX;
    startPanY = e.clientY - panY;
    e.preventDefault();
  }
});

// Wheel / trackpad panning and zooming
editor.addEventListener('wheel', (e) => {
  e.preventDefault();

  if (e.ctrlKey || e.metaKey) {
    // Zoom mode
    const delta = -e.deltaY * 0.001;
    const oldScale = scale;
    scale = Math.max(MIN_SCALE, Math.min(MAX_SCALE, scale + delta));

    // Calculate mouse position relative to editor
    const rect = editor.getBoundingClientRect();
    const mouseX = e.clientX - rect.left;
    const mouseY = e.clientY - rect.top;

    // Adjust pan to keep mouse position fixed
    const scaleDiff = scale / oldScale;
    panX = mouseX - (mouseX - panX) * scaleDiff;
    panY = mouseY - (mouseY - panY) * scaleDiff;
  } else {
    // Pan mode
    panX -= e.deltaX;
    panY -= e.deltaY;
  }

  applyViewportTransform();
  updateConnections();
}, { passive: false });

window.addEventListener('mousemove', (e) => {
  if (isPanning) {
    panX = e.clientX - startPanX;
    panY = e.clientY - startPanY;
    applyViewportTransform();
    updateConnections();
  }
  if (isDraggingNode && draggedNode) {
    const canvasRect = canvas.getBoundingClientRect();
    draggedNode.style.left = `${(e.clientX - canvasRect.left) / scale - nodeOffsetX}px`;
    draggedNode.style.top = `${(e.clientY - canvasRect.top) / scale - nodeOffsetY}px`;
    updateConnections();
  }
  if (isResizing && resizeNode) {
    const dw = (e.clientX - resizeStartX) / scale;
    const dh = (e.clientY - resizeStartY) / scale;
    resizeNode.style.width = Math.max(240, resizeStartW + dw) + 'px';
    resizeNode.style.height = Math.max(100, resizeStartH + dh) + 'px';
    const plotEl = resizeNode.querySelector('.plot-container');
    if (plotEl) Plotly.Plots.resize(plotEl);
    updateConnections();
  }
  if (isWiring && tempWire && wireStartSocket) {
    const canvasRect = canvas.getBoundingClientRect();
    const mx = (e.clientX - canvasRect.left) / scale;
    const my = (e.clientY - canvasRect.top) / scale;
    const sr = getSocketCenter(wireStartSocket);
    if (wireStartIsOutput) {
      tempWire.setAttribute('d', bezierPath(sr.x, sr.y, mx, my));
    } else {
      tempWire.setAttribute('d', bezierPath(mx, my, sr.x, sr.y));
    }
  }
});

window.addEventListener('mouseup', (e) => {
  if (isPanning) isPanning = false;
  if (isDraggingNode) { isDraggingNode = false; draggedNode = null; }
  if (isResizing) { isResizing = false; resizeNode = null; }
  if (isWiring) {
    if (tempWire) { tempWire.remove(); tempWire = null; }
    isWiring = false;
    wireStartSocket = null;
  }
});

editor.addEventListener('contextmenu', (e) => e.preventDefault());

// ===== Node Dragging (via headers) =====
document.addEventListener('mousedown', (e) => {
  const header = e.target.closest('.node-header');
  if (!header || e.button !== 0) return;
  const node = header.closest('.node');
  isDraggingNode = true;
  draggedNode = node;
  const canvasRect = canvas.getBoundingClientRect();
  const nodeLeft = parseFloat(node.style.left || 0);
  const nodeTop = parseFloat(node.style.top || 0);
  nodeOffsetX = (e.clientX - canvasRect.left) / scale - nodeLeft;
  nodeOffsetY = (e.clientY - canvasRect.top) / scale - nodeTop;
  node.style.zIndex = 20;
  document.querySelectorAll('.node').forEach(n => { if (n !== node) n.style.zIndex = 10; });
  e.preventDefault();
});

// ===== Node Resizing =====
document.addEventListener('mousedown', (e) => {
  const handle = e.target.closest('.node-resize');
  if (!handle || e.button !== 0) return;
  const node = handle.closest('.node');
  isResizing = true;
  resizeNode = node;
  resizeStartX = e.clientX;
  resizeStartY = e.clientY;
  resizeStartW = node.offsetWidth;
  resizeStartH = node.offsetHeight;
  e.preventDefault();
  e.stopPropagation();
});

// ===== Socket Wiring =====
document.addEventListener('mousedown', (e) => {
  const socket = e.target.closest('.socket');
  if (!socket || e.button !== 0) return;

  if (socket.classList.contains('output')) {
    // Start wiring from output
    isWiring = true;
    wireStartSocket = socket;
    wireStartIsOutput = true;
    const sr = getSocketCenter(socket);
    tempWire = document.createElementNS('http://www.w3.org/2000/svg', 'path');
    tempWire.classList.add('wire', 'active');
    tempWire.style.stroke = PORT_COLORS[socket.dataset.port] || '#fff';
    tempWire.setAttribute('d', bezierPath(sr.x, sr.y, sr.x, sr.y));
    svgLayer.appendChild(tempWire);
    e.preventDefault();
    e.stopPropagation();
  } else if (socket.classList.contains('input')) {
    const nodeId = socket.dataset.node;
    const port = socket.dataset.port;
    const existing = connections.find(c => c.toNode === nodeId && c.toPort === port);

    if (existing) {
      // Disconnect existing wire and start re-dragging from the output end
      connections = connections.filter(c => c !== existing);
      updateConnections();
      // Start wiring from the original output socket
      const fromSocket = document.querySelector(`#${existing.fromNode} .socket.output[data-port="${existing.fromPort}"]`);
      if (fromSocket) {
        isWiring = true;
        wireStartSocket = fromSocket;
        wireStartIsOutput = true;
        const sr = getSocketCenter(fromSocket);
        tempWire = document.createElementNS('http://www.w3.org/2000/svg', 'path');
        tempWire.classList.add('wire', 'active');
        tempWire.style.stroke = PORT_COLORS[fromSocket.dataset.port] || '#fff';
        tempWire.setAttribute('d', bezierPath(sr.x, sr.y, sr.x, sr.y));
        svgLayer.appendChild(tempWire);
      }
    } else {
      // No existing connection, start wiring from input
      isWiring = true;
      wireStartSocket = socket;
      wireStartIsOutput = false;
      const sr = getSocketCenter(socket);
      tempWire = document.createElementNS('http://www.w3.org/2000/svg', 'path');
      tempWire.classList.add('wire', 'active');
      tempWire.style.stroke = PORT_COLORS[socket.dataset.port] || '#fff';
      tempWire.setAttribute('d', bezierPath(sr.x, sr.y, sr.x, sr.y));
      svgLayer.appendChild(tempWire);
    }
    e.preventDefault();
    e.stopPropagation();
  }
});

document.addEventListener('mouseup', (e) => {
  if (!isWiring || !wireStartSocket) return;
  const socket = e.target.closest('.socket');
  if (socket && socket !== wireStartSocket) {
    let fromSocket, toSocket;
    if (wireStartIsOutput && socket.classList.contains('input')) {
      fromSocket = wireStartSocket;
      toSocket = socket;
    } else if (!wireStartIsOutput && socket.classList.contains('output')) {
      fromSocket = socket;
      toSocket = wireStartSocket;
    }
    if (fromSocket && toSocket) {
      const fromPort = fromSocket.dataset.port;
      const toPort = toSocket.dataset.port;
      // Validate port type compatibility
      if (fromPort === toPort) {
        const fromNode = fromSocket.dataset.node;
        const toNode = toSocket.dataset.node;
        // No self-connections
        if (fromNode !== toNode) {
          // Remove existing connection to this input (one input = one wire)
          connections = connections.filter(c => !(c.toNode === toNode && c.toPort === toPort));
          connections.push({ fromNode, fromPort, toNode, toPort });
          updateConnections();

          // Auto-populate config nodes when connected
          const toNodeInfo = nodeRegistry[toNode];
          if (toNodeInfo && toNodeInfo.type) {
            const typeDef = NODE_TYPES[toNodeInfo.type];
            if (typeDef && typeDef.category === 'parameter' && typeDef.execute) {
              // Execute the config node to populate dropdowns/options
              // Check if we have the necessary data before executing
              const shouldExecute =
                (toPort === 'model' && state.sessionId) || // Has model data
                (toPort === 'reactions'); // Has reactions data

              if (shouldExecute) {
                setTimeout(() => {
                  typeDef.execute(toNode).catch(e => {
                    console.error(`Failed to auto-populate ${toNode}:`, e);
                  });
                }, 100);
              }
            }
          }
        }
      } else {
        showToast(`Port mismatch: ${fromPort} ≠ ${toPort}`);
      }
    }
  }
});

// ===== Connection Drawing =====
function getSocketCenter(socket) {
  const rect = socket.getBoundingClientRect();
  const canvasRect = canvas.getBoundingClientRect();
  return {
    x: (rect.left + rect.width / 2 - canvasRect.left) / scale,
    y: (rect.top + rect.height / 2 - canvasRect.top) / scale,
  };
}

function bezierPath(x1, y1, x2, y2) {
  const dx = Math.abs(x2 - x1) * 0.5;
  return `M ${x1} ${y1} C ${x1 + dx} ${y1}, ${x2 - dx} ${y2}, ${x2} ${y2}`;
}

function updateConnections() {
  // Store transmitting state before removing wires
  const transmittingWires = new Set();
  svgLayer.querySelectorAll('.wire.connected.transmitting').forEach(w => {
    const id = w.getAttribute('id');
    if (id) transmittingWires.add(id);
  });

  svgLayer.querySelectorAll('.wire.connected').forEach(w => w.remove());
  // Reset all socket connected state
  document.querySelectorAll('.socket.connected').forEach(s => s.classList.remove('connected'));
  connections.forEach(conn => {
    const fromSocket = document.querySelector(`#${conn.fromNode} .socket.output[data-port="${conn.fromPort}"]`);
    const toSocket = document.querySelector(`#${conn.toNode} .socket.input[data-port="${conn.toPort}"]`);
    if (!fromSocket || !toSocket) return;
    fromSocket.classList.add('connected');
    toSocket.classList.add('connected');
    const from = getSocketCenter(fromSocket);
    const to = getSocketCenter(toSocket);
    const path = document.createElementNS('http://www.w3.org/2000/svg', 'path');
    const wireId = `wire-${conn.fromNode}-${conn.toNode}`;
    path.classList.add('wire', 'connected');
    path.setAttribute('id', wireId);
    path.setAttribute('d', bezierPath(from.x, from.y, to.x, to.y));
    path.setAttribute('data-port-type', conn.fromPort);
    path.style.stroke = PORT_COLORS[conn.fromPort] || '#888';

    // Restore transmitting state if it was active
    if (transmittingWires.has(wireId)) {
      path.classList.add('transmitting');
    }

    svgLayer.appendChild(path);
  });
}

// ===== Node Factory =====
function createNode(nodeType, x, y) {
  const typeDef = NODE_TYPES[nodeType];
  if (!typeDef) { console.error('Unknown node type:', nodeType); return null; }

  nodeIdCounter++;
  const nodeId = `node-${nodeIdCounter}`;

  const node = document.createElement('div');
  const isLargeNode = ['viewer', 'result', 'parameter'].includes(typeDef.category);
  node.className = `node${isLargeNode ? ' viewer' : ''}`;
  node.id = nodeId;
  node.dataset.type = typeDef.category;
  node.dataset.nodeType = nodeType;
  node.style.left = `${x}px`;
  node.style.top = `${y}px`;
  if (typeDef.defaultWidth) node.style.width = `${typeDef.defaultWidth}px`;

  // Header
  const header = document.createElement('div');
  header.className = `node-header ${typeDef.headerClass}`;
  header.innerHTML = `
    <span>${typeDef.title}</span>
    <button class="btn-close" onclick="removeNode('${nodeId}')">&times;</button>
  `;
  node.appendChild(header);

  // Body
  const body = document.createElement('div');
  body.className = 'node-body';

  // Input sockets
  typeDef.inputs.forEach(inp => {
    body.innerHTML += `
      <div class="socket-row left">
        <div class="socket input" data-node="${nodeId}" data-port="${inp.port}"></div>
        <span class="socket-label">${inp.label}</span>
      </div>
    `;
  });

  // Custom body content
  if (typeDef.createBody) {
    body.innerHTML += typeDef.createBody(nodeId);
  }

  // Output sockets
  typeDef.outputs.forEach(out => {
    body.innerHTML += `
      <div class="socket-row right">
        <span class="socket-label">${out.label}</span>
        <div class="socket output" data-node="${nodeId}" data-port="${out.port}"></div>
      </div>
    `;
  });

  node.appendChild(body);

  // Resize handle
  const resize = document.createElement('div');
  resize.className = 'node-resize';
  node.appendChild(resize);

  canvas.appendChild(node);

  nodeRegistry[nodeId] = { type: nodeType, el: node, data: {} };

  // Run init hook
  if (typeDef.onInit) typeDef.onInit(nodeId);

  return nodeId;
}

function removeNode(nodeId) {
  const el = document.getElementById(nodeId);
  if (el) el.remove();
  connections = connections.filter(c => c.fromNode !== nodeId && c.toNode !== nodeId);
  delete nodeRegistry[nodeId];
  cleanupPlotResize(nodeId);
  updateConnections();
}

// ===== Node Loading State =====
function setNodeLoading(nodeId, loading) {
  const el = document.getElementById(nodeId);
  if (!el) return;
  if (loading) {
    el.classList.add('loading');
    // Mark all input wires as transmitting
    const inputConns = connections.filter(c => c.toNode === nodeId);
    inputConns.forEach(conn => {
      const wireId = `wire-${conn.fromNode}-${conn.toNode}`;
      const wire = document.getElementById(wireId);
      if (wire) wire.classList.add('transmitting');
    });
  } else {
    el.classList.remove('loading');
    // Remove transmitting state from all input wires
    const inputConns = connections.filter(c => c.toNode === nodeId);
    inputConns.forEach(conn => {
      const wireId = `wire-${conn.fromNode}-${conn.toNode}`;
      const wire = document.getElementById(wireId);
      if (wire) wire.classList.remove('transmitting');
    });
  }
}

// ===== Auto-Chain Generation =====

// Find an existing chain ending with a model-builder that has a model output
function findExistingModelBuilder() {
  for (const [id, info] of Object.entries(nodeRegistry)) {
    if (info.type === 'model-builder') {
      // Check if this model-builder is connected to a reaction-network
      const conn = connections.find(c => c.toNode === id && c.toPort === 'reactions');
      if (conn && nodeRegistry[conn.fromNode]?.type === 'reaction-network') {
        return { modelBuilderId: id, reactionNetworkId: conn.fromNode };
      }
    }
  }
  return null;
}

function findExistingReactionNetwork() {
  for (const [id, info] of Object.entries(nodeRegistry)) {
    if (info.type === 'reaction-network') return id;
  }
  return null;
}

function getNodePosition(nodeId) {
  const el = document.getElementById(nodeId);
  if (!el) return { x: 100, y: 150 };
  return { x: parseFloat(el.style.left) || 0, y: parseFloat(el.style.top) || 0 };
}

function getNodeSize(nodeId) {
  const el = document.getElementById(nodeId);
  if (!el) return { w: 260, h: 200 };
  return { w: el.offsetWidth, h: el.offsetHeight };
}

// Count how many viewers are already attached to a model-builder
function countDownstreamViewers(modelBuilderId) {
  return connections.filter(c => c.fromNode === modelBuilderId && c.fromPort === 'model').length;
}

// Simple collision detection — shift node down if overlapping
function resolveOverlap(x, y, width, height, excludeNodeId) {
  let maxAttempts = 20;
  let curY = y;
  while (maxAttempts-- > 0) {
    let overlaps = false;
    for (const [id, info] of Object.entries(nodeRegistry)) {
      if (id === excludeNodeId) continue;
      const pos = getNodePosition(id);
      const size = getNodeSize(id);
      if (x < pos.x + size.w && x + width > pos.x &&
          curY < pos.y + size.h && curY + height > pos.y) {
        curY = pos.y + size.h + 30;
        overlaps = true;
        break;
      }
    }
    if (!overlaps) break;
  }
  return curY;
}

function addNodeFromMenu(nodeType) {
  closeDropdown();

  // Simple strategy: just create the node at a reasonable position
  const typeDef = NODE_TYPES[nodeType];
  if (!typeDef) return;

  // Find a good position based on existing nodes
  let x = 80;
  let y = 150;

  // If there are existing nodes, place new node to the right
  const existingNodes = Object.keys(nodeRegistry);
  if (existingNodes.length > 0) {
    let maxX = 0;
    for (const id of existingNodes) {
      const pos = getNodePosition(id);
      const size = getNodeSize(id);
      if (pos.x + size.w > maxX) {
        maxX = pos.x + size.w;
      }
    }
    x = maxX + 60;
  }

  const width = typeDef.defaultWidth || 280;
  y = resolveOverlap(x, y, width, 300, null);

  createNode(nodeType, x, y);
}

function addResultNode(nodeType) {
  // This function is no longer used - kept for compatibility
  // All nodes are now added via addNodeFromMenu
  addNodeFromMenu(nodeType);
}

// ===== Quick Add Chain Generation =====
function addQuickAddChain(chainType) {
  closeDropdown();

  // Map legacy node types to their new chain equivalents
  const chainMap = {
    'siso-analysis': { params: 'siso-params', result: 'siso-result' },
    'rop-cloud': { params: 'rop-cloud-params', result: 'rop-cloud-result' },
    'fret-heatmap': { params: 'fret-params', result: 'fret-result' },
    'parameter-scan-1d': { params: 'scan-1d-params', result: 'scan-1d-result' },
    'parameter-scan-2d': { params: 'scan-2d-params', result: 'scan-2d-result' },
    'rop-polyhedron': { params: 'rop-poly-params', result: 'rop-poly-result' },
  };

  const chain = chainMap[chainType];
  if (!chain) {
    console.error('Unknown quick add chain type:', chainType);
    return;
  }

  // Check for existing nodes and reuse them
  let rnId = findExistingReactionNetwork();
  let mbId = null;

  if (!rnId) {
    // No reaction network exists, create one
    rnId = createNode('reaction-network', 80, 150);
  }

  // Check for existing model-builder connected to this reaction network
  const existing = findExistingModelBuilder();
  if (existing && existing.reactionNetworkId === rnId) {
    mbId = existing.modelBuilderId;
  } else {
    // Create model-builder and connect to reaction network
    const rnPos = getNodePosition(rnId);
    const rnSize = getNodeSize(rnId);
    const mbX = rnPos.x + rnSize.w + 60;
    const mbY = resolveOverlap(mbX, rnPos.y, 260, 200, null);
    mbId = createNode('model-builder', mbX, mbY);
    connections.push({ fromNode: rnId, fromPort: 'reactions', toNode: mbId, toPort: 'reactions' });
  }

  // Create params and result nodes
  const mbPos = getNodePosition(mbId);
  const mbSize = getNodeSize(mbId);
  const paramsX = mbPos.x + mbSize.w + 60;
  const nDownstream = countDownstreamViewers(mbId);
  const paramsY = resolveOverlap(paramsX, mbPos.y + nDownstream * 50, 320, 300, null);
  const paramsId = createNode(chain.params, paramsX, paramsY);

  const paramsSize = getNodeSize(paramsId);
  const resultX = paramsX + paramsSize.w + 60;
  const resultY = resolveOverlap(resultX, paramsY, 420, 300, null);
  const resultId = createNode(chain.result, resultX, resultY);

  // Connect them
  connections.push({ fromNode: mbId, fromPort: 'model', toNode: paramsId, toPort: 'model' });
  connections.push({ fromNode: paramsId, fromPort: 'params', toNode: resultId, toPort: 'params' });

  // Special case: ROP cloud params also needs reactions connection
  if (chain.params === 'rop-cloud-params') {
    connections.push({ fromNode: rnId, fromPort: 'reactions', toNode: paramsId, toPort: 'reactions' });
  }

  updateConnections();

  // Auto-populate the params node if model data is available
  const paramsTypeDef = NODE_TYPES[chain.params];
  if (paramsTypeDef && paramsTypeDef.execute) {
    // Check if we have model data or reactions data
    const hasModelData = state.sessionId;
    const hasReactionsData = chain.params === 'rop-cloud-params'; // ROP cloud uses reactions

    if (hasModelData || hasReactionsData) {
      setTimeout(() => {
        paramsTypeDef.execute(paramsId).catch(e => {
          console.error(`Failed to auto-populate ${paramsId}:`, e);
        });
      }, 100);
    }
  }
}

// ===== Toolbar / Dropdown =====
const addNodeBtn = document.getElementById('add-node-btn');
const addNodeMenu = document.getElementById('add-node-menu');
const legacyNodesBtn = document.getElementById('legacy-nodes-btn');
const legacyNodesMenu = document.getElementById('legacy-nodes-menu');

addNodeBtn.addEventListener('click', (e) => {
  e.stopPropagation();
  addNodeMenu.classList.toggle('open');
  legacyNodesMenu.classList.remove('open');
});

legacyNodesBtn.addEventListener('click', (e) => {
  e.stopPropagation();
  legacyNodesMenu.classList.toggle('open');
  addNodeMenu.classList.remove('open');
});

document.addEventListener('click', (e) => {
  if (!addNodeMenu.contains(e.target) && e.target !== addNodeBtn) {
    addNodeMenu.classList.remove('open');
  }
  if (!legacyNodesMenu.contains(e.target) && e.target !== legacyNodesBtn) {
    legacyNodesMenu.classList.remove('open');
  }
});

addNodeMenu.querySelectorAll('.menu-item').forEach(item => {
  item.addEventListener('click', () => {
    addNodeFromMenu(item.dataset.type);
  });
});

legacyNodesMenu.querySelectorAll('.menu-item').forEach(item => {
  item.addEventListener('click', () => {
    addQuickAddChain(item.dataset.type);
  });
});

function closeDropdown() {
  addNodeMenu.classList.remove('open');
  legacyNodesMenu.classList.remove('open');
}

// ===== Toast Notifications =====
function showToast(message, duration = 2500) {
  const container = document.getElementById('toast-container');
  const toast = document.createElement('div');
  toast.className = 'toast';
  toast.textContent = message;
  container.appendChild(toast);
  requestAnimationFrame(() => toast.classList.add('show'));
  setTimeout(() => {
    toast.classList.remove('show');
    setTimeout(() => toast.remove(), 300);
  }, duration);
}

// ===== API Helpers =====
async function api(endpoint, data) {
  setStatus('working', 'Computing...');
  try {
    const resp = await fetch(`${API}/api/${endpoint}`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(data),
    });

    // Check if response is JSON
    const contentType = resp.headers.get('content-type');
    if (!contentType || !contentType.includes('application/json')) {
      throw new Error('Backend server not responding. Please ensure Julia server is running.');
    }

    const json = await resp.json();
    if (json.error) throw new Error(json.error);
    setStatus('done', 'Done');
    return json;
  } catch (e) {
    setStatus('error', e.message);
    throw e;
  }
}

function setStatus(cls, text) {
  const badge = document.getElementById('status-badge');
  badge.className = `badge ${cls}`;
  badge.textContent = text;
  if (cls === 'done') setTimeout(() => {
    badge.className = 'badge';
    badge.textContent = 'Ready';
  }, 3000);
}

// ===== Reaction Editor =====
function getReactionsFromNode(nodeId) {
  const list = document.getElementById(`${nodeId}-reactions-list`);
  if (!list) return { reactions: [], kds: [] };
  const rows = list.querySelectorAll('.reaction-row');
  const reactions = [];
  const kds = [];
  rows.forEach(row => {
    const rule = row.querySelector('.reaction-input').value.trim();
    const kd = parseFloat(row.querySelector('.kd-input').value);
    if (rule) {
      reactions.push(rule);
      kds.push(Number.isFinite(kd) ? kd : null);
    }
  });
  return { reactions, kds };
}

function addReactionRow(nodeId, rule = '', kd = 1e-3) {
  const list = document.getElementById(`${nodeId}-reactions-list`);
  if (!list) return;
  const row = document.createElement('div');
  row.className = 'reaction-row';
  row.innerHTML = `
    <input type="text" class="reaction-input" value="${rule}" placeholder="A + B <-> C">
    <input type="number" class="kd-input" value="${kd == null ? '' : kd}" step="any" min="1e-12" placeholder="optional">
    <button class="btn-remove" title="Remove">&times;</button>
  `;

  const removeBtn = row.querySelector('.btn-remove');
  removeBtn.onclick = () => {
    row.remove();
    triggerAutoModelBuild(nodeId);
  };

  // Add event listeners for auto-build
  const reactionInput = row.querySelector('.reaction-input');
  const kdInput = row.querySelector('.kd-input');

  [reactionInput, kdInput].forEach(input => {
    input.addEventListener('input', () => {
      clearTimeout(input._autoTimer);
      input._autoTimer = setTimeout(() => {
        triggerAutoModelBuild(nodeId);
      }, 1000);
    });
  });

  list.appendChild(row);
}

// ===== Build Model =====
async function buildModel(modelBuilderNodeId) {
  // Find connected reaction-network
  const conn = connections.find(c => c.toNode === modelBuilderNodeId && c.toPort === 'reactions');
  if (!conn) {
    showToast('Model Builder has no Reaction Network connected');
    return;
  }
  const rnNodeId = conn.fromNode;
  const { reactions, kds } = getReactionsFromNode(rnNodeId);
  if (reactions.length === 0) { showToast('Add at least one reaction'); return; }
  if (kds.some(kd => kd == null || kd <= 0)) { showToast('Model Builder requires Kd for every reaction (> 0)'); return; }

  setNodeLoading(modelBuilderNodeId, true);
  try {
    const data = await api('build_model', { reactions, kd: kds });
    state.sessionId = data.session_id;
    state.model = data;
    state.qK_syms = [...data.q_sym, ...data.K_sym];

    // Update model info display
    const infoEl = document.getElementById(`${modelBuilderNodeId}-model-info`);
    const infoText = document.getElementById(`${modelBuilderNodeId}-model-info-text`);
    if (infoEl && infoText) {
      const info = `n=${data.n}, d=${data.d}, r=${data.r}\nSpecies: ${data.x_sym.join(', ')}\nTotals: ${data.q_sym.join(', ')}\nConstants: ${data.K_sym.join(', ')}`;
      infoEl.style.display = '';
      infoText.textContent = info;
    }

    // Store model builder node reference
    nodeRegistry[modelBuilderNodeId].data.built = true;

    showToast('Model built successfully');

    // Trigger all downstream viewers
    onModelBuilt(modelBuilderNodeId);
  } catch (e) {
    console.error('Build model failed:', e);
  }
  setNodeLoading(modelBuilderNodeId, false);
}

// ===== Downstream Viewer Auto-Execution =====
async function onModelBuilt(modelBuilderNodeId) {
  const downstream = connections.filter(c => c.fromNode === modelBuilderNodeId && c.fromPort === 'model');

  for (const conn of downstream) {
    const viewerInfo = nodeRegistry[conn.toNode];
    if (!viewerInfo) continue;
    const typeDef = NODE_TYPES[viewerInfo.type];
    if (typeDef && typeDef.execute) {
      // Mark input wires as transmitting
      const wireId = `wire-${conn.fromNode}-${conn.toNode}`;
      const wire = document.getElementById(wireId);
      if (wire) wire.classList.add('transmitting');

      // Execute asynchronously (don't await — run in parallel)
      typeDef.execute(conn.toNode).catch(e => {
        console.error(`Node ${conn.toNode} (${viewerInfo.type}) failed:`, e);
      }).finally(() => {
        // Remove transmitting state after a short delay
        setTimeout(() => {
          if (wire) wire.classList.remove('transmitting');
        }, 500);
      });
    }
  }
}

// ===== Recompute Functions =====
function recomputeSISO(nodeId) {
  const typeDef = NODE_TYPES['siso-analysis'];
  if (typeDef.execute) typeDef.execute(nodeId);
}

async function computeSISOResult(nodeId) {
  // Find the connected params node
  const paramsConn = connections.find(c => c.toNode === nodeId && c.toPort === 'params');
  if (!paramsConn) {
    showToast('Connect a SISO Config node first');
    return;
  }

  const paramsNode = nodeRegistry[paramsConn.fromNode];
  if (!paramsNode || !paramsNode.data || !paramsNode.data.config) {
    showToast('SISO Config node has no configuration');
    return;
  }

  const config = paramsNode.data.config;
  setNodeLoading(nodeId, true);
  const contentEl = document.getElementById(`${nodeId}-content`);

  try {
    const data = await api('siso_paths', {
      session_id: state.sessionId,
      change_qK: config.changeQK || config.change_qK
    });

    let html = `<div style="margin-bottom:8px;"><strong>${data.n_paths}</strong> paths, <strong>${data.sources.length}</strong> sources, <strong>${data.sinks.length}</strong> sinks</div>`;
    html += '<div class="path-list">';
    data.paths.forEach(p => {
      const permStr = p.perms.map(pr => `[${pr.join(',')}]`).join(' → ');
      html += `<div class="path-item" data-idx="${p.idx}" data-qk="${config.changeQK || config.change_qK}" data-node="${nodeId}" onclick="selectSISOPath(this)">#${p.idx}: ${permStr}</div>`;
    });
    html += '</div>';
    html += `<div class="plot-container" id="${nodeId}-traj-plot" style="display:none;"></div>`;
    contentEl.innerHTML = html;
  } catch (e) {
    contentEl.innerHTML = `<div class="node-error">${e.message}</div>`;
  } finally {
    setNodeLoading(nodeId, false);
  }
}

function recomputeROPCloud(nodeId) {
  const typeDef = NODE_TYPES['rop-cloud'];
  if (typeDef.execute) typeDef.execute(nodeId);
}

function recomputeHeatmap(nodeId) {
  const typeDef = NODE_TYPES['fret-heatmap'];
  if (typeDef.execute) typeDef.execute(nodeId);
}

function parseSpeciesFromReactionSide(side) {
  const species = [];
  side.split('+').forEach(term => {
    const t = term.trim();
    if (!t) return;
    const m = t.match(/^([0-9]+)?\s*([A-Za-z][A-Za-z0-9_]*)$/);
    if (m) species.push(m[2]);
  });
  return species;
}

function inferSpeciesOrderFromReactions(reactions) {
  const allSet = new Set();
  const productSet = new Set();
  reactions.forEach(rule => {
    const m = rule.match(/<->|<=>|↔/);
    if (!m) return;
    const parts = rule.split(m[0]);
    if (parts.length !== 2) return;
    const left = parseSpeciesFromReactionSide(parts[0]);
    const right = parseSpeciesFromReactionSide(parts[1]);
    left.forEach(s => allSet.add(s));
    right.forEach(s => {
      allSet.add(s);
      productSet.add(s);
    });
  });

  const allSpecies = Array.from(allSet).sort();
  const productSpecies = Array.from(productSet).sort();
  const freeSpecies = allSpecies.filter(s => !productSet.has(s));
  const orderedSpecies = [...freeSpecies, ...productSpecies];
  return { species: orderedSpecies, productSpecies };
}

function refreshROPCloudTargetOptions(nodeId, reactions = null) {
  const sel = document.getElementById(`${nodeId}-target-species`);
  if (!sel) return;

  if (!reactions) {
    const rxConn = connections.find(c => c.toNode === nodeId && c.toPort === 'reactions');
    if (rxConn) reactions = getReactionsFromNode(rxConn.fromNode).reactions;
  }
  reactions = reactions || [];

  const { species, productSpecies } = inferSpeciesOrderFromReactions(reactions);
  const preferred = productSpecies.length ? productSpecies : species;
  const orderedTargets = [...preferred, ...species.filter(s => !preferred.includes(s))];

  const prev = sel.value;
  sel.innerHTML = '';
  if (!orderedTargets.length) {
    const opt = document.createElement('option');
    opt.value = '';
    opt.textContent = '(target)';
    sel.appendChild(opt);
    return;
  }

  orderedTargets.forEach(sym => {
    const opt = document.createElement('option');
    opt.value = sym;
    opt.textContent = sym;
    sel.appendChild(opt);
  });

  if (prev && orderedTargets.includes(prev)) {
    sel.value = prev;
  }
}

function updateROPCloudMode(nodeId) {
  const mode = document.getElementById(`${nodeId}-sampling-mode`)?.value || 'x_space';
  const xParams = document.getElementById(`${nodeId}-xspace-params`);
  const qkParams = document.getElementById(`${nodeId}-qk-params`);
  if (xParams) xParams.style.display = mode === 'x_space' ? '' : 'none';
  if (qkParams) qkParams.style.display = mode === 'qk' ? '' : 'none';
  if (mode === 'x_space') refreshROPCloudTargetOptions(nodeId);
}

// ===== SISO Path Selection =====
async function selectSISOPath(el) {
  // Deselect siblings
  const parent = el.closest('.path-list');
  if (parent) parent.querySelectorAll('.path-item').forEach(p => p.classList.remove('selected'));
  el.classList.add('selected');

  const pathIdx = parseInt(el.dataset.idx);
  const changeQK = el.dataset.qk;
  const nodeId = el.dataset.node;

  try {
    const data = await api('siso_trajectory', {
      session_id: state.sessionId,
      change_qK: changeQK,
      path_idx: pathIdx,
    });
    const plotEl = document.getElementById(`${nodeId}-traj-plot`);
    if (plotEl) {
      plotEl.style.display = '';
      plotTrajectory(data, `${nodeId}-traj-plot`);
    }
  } catch (e) {
    console.error('Trajectory failed:', e);
  }
}

// ===== Plotly Renderers =====

function plotRegimeGraph(data, plotId) {
  const { nodes, edges } = data;
  const n = nodes.length;

  const positions = {};
  nodes.forEach((node, i) => {
    const angle = (2 * Math.PI * i) / n;
    const r = 2 + Math.sqrt(n) * 0.5;
    positions[node.id] = { x: r * Math.cos(angle), y: r * Math.sin(angle) };
  });

  // Spring layout
  for (let iter = 0; iter < 200; iter++) {
    const forces = {};
    nodes.forEach(nd => { forces[nd.id] = { x: 0, y: 0 }; });
    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) {
        const a = nodes[i].id, b = nodes[j].id;
        let dx = positions[a].x - positions[b].x;
        let dy = positions[a].y - positions[b].y;
        let dist = Math.sqrt(dx * dx + dy * dy) + 0.01;
        let f = 3.0 / (dist * dist);
        forces[a].x += f * dx / dist; forces[a].y += f * dy / dist;
        forces[b].x -= f * dx / dist; forces[b].y -= f * dy / dist;
      }
    }
    edges.forEach(e => {
      let dx = positions[e.target].x - positions[e.source].x;
      let dy = positions[e.target].y - positions[e.source].y;
      let dist = Math.sqrt(dx * dx + dy * dy) + 0.01;
      let f = dist * 0.1;
      forces[e.source].x += f * dx / dist; forces[e.source].y += f * dy / dist;
      forces[e.target].x -= f * dx / dist; forces[e.target].y -= f * dy / dist;
    });
    const cooling = 0.95 - 0.5 * (iter / 200);
    nodes.forEach(nd => {
      positions[nd.id].x += forces[nd.id].x * cooling;
      positions[nd.id].y += forces[nd.id].y * cooling;
    });
  }

  const traces = [];

  // Add edges with hover info
  edges.forEach(e => {
    const sourceNode = nodes.find(n => n.id === e.source);
    const targetNode = nodes.find(n => n.id === e.target);
    traces.push({
      x: [positions[e.source].x, positions[e.target].x],
      y: [positions[e.source].y, positions[e.target].y],
      mode: 'lines',
      line: { color: '#444', width: 1 },
      hoverinfo: 'text',
      text: `Edge: #${e.source} → #${e.target}<br>From: [${sourceNode.perm.join(',')}]<br>To: [${targetNode.perm.join(',')}]`,
      showlegend: false,
      type: 'scatter',
    });
  });

  const colors = nodes.map(nd => {
    if (nd.singular) return '#ff6b6b';
    if (!nd.asymptotic) return '#ffd43b';
    return '#51cf66';
  });

  traces.push({
    x: nodes.map(nd => positions[nd.id].x),
    y: nodes.map(nd => positions[nd.id].y),
    mode: 'markers+text', type: 'scatter',
    marker: { size: 18, color: colors, line: { width: 1, color: '#666' } },
    text: nodes.map(nd => `${nd.id}`),
    textposition: 'middle center',
    textfont: { size: 9, color: '#fff' },
    hovertext: nodes.map(nd =>
      `#${nd.id} [${nd.perm.join(',')}]<br>${nd.asymptotic ? 'Asymptotic' : 'Non-Asymp'} ${nd.singular ? 'Singular' : 'Invertible'}<br>nullity=${nd.nullity}`
    ),
    hoverinfo: 'text',
  });

  const layout = {
    autosize: true,
    showlegend: false,
    xaxis: { visible: false }, yaxis: { visible: false, scaleanchor: 'x' },
    paper_bgcolor: '#111', plot_bgcolor: '#111',
    margin: { t: 40, b: 20, l: 20, r: 20 },
    title: { text: `${nodes.length} vertices, ${edges.length} edges`, font: { color: '#888', size: 11 }, y: 0.98, yanchor: 'top' },
  };

  Plotly.newPlot(plotId, traces, layout, { responsive: true, displayModeBar: false });
}

function plotTrajectory(data, plotId) {
  const { change_values, logx, regimes, x_sym, change_sym } = data;
  const nSpecies = x_sym.length;
  const nPoints = change_values.length;

  const uniqueRegimes = [...new Set(regimes)];
  const palette = ['#6c8cff', '#51cf66', '#ff6b6b', '#ffd43b', '#4ecdc4', '#e599f7', '#ff922b', '#74c0fc', '#f06595', '#a9e34b'];
  const regimeColor = {};
  uniqueRegimes.forEach((r, i) => { regimeColor[r] = palette[i % palette.length]; });

  const traces = [];
  for (let s = 0; s < nSpecies; s++) {
    let i = 0;
    while (i < nPoints) {
      const regime = regimes[i];
      const segX = [], segY = [];
      while (i < nPoints && regimes[i] === regime) {
        segX.push(change_values[i]);
        segY.push(logx[i][s]);
        i++;
      }
      traces.push({
        x: segX, y: segY, mode: 'lines', type: 'scatter',
        line: { color: regimeColor[regime], width: 2 },
        name: `${x_sym[s]} (rgm ${regime})`,
        legendgroup: x_sym[s],
        showlegend: segX[0] === change_values[0],
      });
    }
  }

  const layout = {
    showlegend: true,
    paper_bgcolor: '#111', plot_bgcolor: '#1a1a1a',
    font: { color: '#888', size: 10 },
    margin: { t: 40, b: 60, l: 70, r: 20 },
    title: { text: `Changing ${change_sym}`, font: { color: '#888', size: 11 }, y: 0.98, yanchor: 'top' },
    xaxis: { title: `log ${change_sym}`, gridcolor: '#333', zerolinecolor: '#444' },
    yaxis: { title: 'log(x)', gridcolor: '#333', zerolinecolor: '#444' },
    legend: { bgcolor: 'rgba(0,0,0,0)', font: { color: '#888', size: 9 } },
  };

  Plotly.newPlot(plotId, traces, layout, { responsive: true, displayModeBar: false });
}

function plotROPCloud(data, plotId) {
  const { reaction_orders, fret_values, q_sym, d } = data;
  const darkLayout = {
    autosize: true,
    paper_bgcolor: '#111', plot_bgcolor: '#1a1a1a',
    font: { color: '#888', size: 10 },
    margin: { t: 40, b: 60, l: 70, r: 20 },
  };

  if (d === 2) {
    const x = reaction_orders.map(r => r[0]);
    const y = reaction_orders.map(r => r[1]);
    const traces = [{
      x, y, mode: 'markers', type: 'scatter',
      marker: {
        size: 2, color: fret_values.map(v => Math.log10(v + 1e-30)),
        colorscale: 'Viridis', showscale: true,
        colorbar: { title: 'log(FRET)', titlefont: { color: '#888', size: 9 }, tickfont: { size: 8 } },
      },
    }];
    Plotly.newPlot(plotId, traces, {
      ...darkLayout,
      title: { text: 'ROP Cloud (2D)', font: { color: '#888', size: 11 }, y: 0.98, yanchor: 'top' },
      xaxis: { title: `\u2202log/\u2202log ${q_sym[0]}`, gridcolor: '#333' },
      yaxis: { title: `\u2202log/\u2202log ${q_sym[1]}`, gridcolor: '#333' },
    }, { responsive: true, displayModeBar: false });
  } else if (d === 3) {
    const x = reaction_orders.map(r => r[0]);
    const y = reaction_orders.map(r => r[1]);
    const z = reaction_orders.map(r => r[2]);
    const traces = [{
      x, y, z, mode: 'markers', type: 'scatter3d',
      marker: {
        size: 2, color: fret_values.map(v => Math.log10(v + 1e-30)),
        colorscale: 'Viridis', showscale: true,
      },
    }];
    Plotly.newPlot(plotId, traces, {
      ...darkLayout,
        title: { text: 'ROP Cloud (3D)', font: { color: '#888', size: 11 }, y: 0.98, yanchor: 'top' },
      scene: {
        xaxis: { title: `\u2202log/\u2202log ${q_sym[0]}`, gridcolor: '#333' },
        yaxis: { title: `\u2202log/\u2202log ${q_sym[1]}`, gridcolor: '#333' },
        zaxis: { title: `\u2202log/\u2202log ${q_sym[2]}`, gridcolor: '#333' },
        bgcolor: '#1a1a1a',
      },
    }, { responsive: true, displayModeBar: false });
  } else {
    const x = reaction_orders.map(r => r[0]);
    const y = reaction_orders.map(r => r[1]);
    const traces = [{
      x, y, mode: 'markers', type: 'scatter',
      marker: { size: 2, color: '#6c8cff', opacity: 0.3 },
    }];
    Plotly.newPlot(plotId, traces, {
      ...darkLayout,
      title: { text: `ROP Cloud (first 2 of ${d} dims)`, font: { color: '#888', size: 11 }, y: 0.98, yanchor: 'top' },
      xaxis: { title: `\u2202log/\u2202log ${q_sym[0]}`, gridcolor: '#333' },
      yaxis: { title: `\u2202log/\u2202log ${q_sym[1]}`, gridcolor: '#333' },
    }, { responsive: true, displayModeBar: false });
  }
}

function plotHeatmap(data, plotId) {
  const { logq1, logq2, fret, regime, bounds, q_sym } = data;
  const logFret = fret.map(row => row.map(v => Math.log10(v + 1e-30)));

  const traces = [
    {
      z: logFret, x: logq1, y: logq2, type: 'heatmap', colorscale: 'Viridis',
      colorbar: { title: 'log(FRET)', titlefont: { color: '#888', size: 9 }, tickfont: { size: 8 } },
    },
    {
      z: bounds, x: logq1, y: logq2, type: 'contour',
      contours: { start: 0.5, end: 0.5, size: 1, coloring: 'none' },
      line: { color: '#fff', width: 2 }, showscale: false,
    },
  ];

  const layout = {
    autosize: true,
    paper_bgcolor: '#111', plot_bgcolor: '#1a1a1a',
    font: { color: '#888', size: 10 },
    margin: { t: 40, b: 60, l: 70, r: 20 },
    title: { text: 'FRET + Regime Boundaries', font: { color: '#888', size: 11 }, y: 0.98, yanchor: 'top' },
    xaxis: { title: `log(${q_sym[0]})`, gridcolor: '#333' },
    yaxis: { title: `log(${q_sym[1]})`, gridcolor: '#333' },
  };

  Plotly.newPlot(plotId, traces, layout, { responsive: true, displayModeBar: false });
}

// ===== Reset View =====
function resetView() {
  panX = 0; panY = 0; scale = 1.0;
  applyViewportTransform();
  updateConnections();
}

// ===== Save / Load Workspace =====

function getNodeSerialData(nodeId, type) {
  switch (type) {
    case 'reaction-network': {
      const { reactions, kds } = getReactionsFromNode(nodeId);
      return { reactions: reactions.map((rule, i) => ({ rule, kd: kds[i] })) };
    }
    case 'siso-params': {
      const changeQK = document.getElementById(`${nodeId}-siso-select`)?.value || '';
      const min = parseFloat(document.getElementById(`${nodeId}-min`)?.value || '-6');
      const max = parseFloat(document.getElementById(`${nodeId}-max`)?.value || '6');
      return { changeQK, min, max, config: nodeRegistry[nodeId]?.data?.config };
    }
    case 'rop-cloud': {
      const mode = document.getElementById(`${nodeId}-sampling-mode`)?.value || 'x_space';
      const samples = parseInt(document.getElementById(`${nodeId}-samples`)?.value || '10000');
      const span = parseInt(document.getElementById(`${nodeId}-span`)?.value || '6');
      const logxMin = parseFloat(document.getElementById(`${nodeId}-logx-min`)?.value || '-6');
      const logxMax = parseFloat(document.getElementById(`${nodeId}-logx-max`)?.value || '6');
      const targetSpecies = document.getElementById(`${nodeId}-target-species`)?.value || '';
      return { mode, samples, span, logxMin, logxMax, targetSpecies };
    }
    case 'fret-heatmap': {
      const grid = parseInt(document.getElementById(`${nodeId}-grid`)?.value || '80');
      return { grid };
    }
    case 'scan-1d-params': {
      const param = document.getElementById(`${nodeId}-param`)?.value || '';
      const min = parseFloat(document.getElementById(`${nodeId}-min`)?.value || '-6');
      const max = parseFloat(document.getElementById(`${nodeId}-max`)?.value || '6');
      const points = parseInt(document.getElementById(`${nodeId}-points`)?.value || '200');
      const expr = (document.getElementById(`${nodeId}-expr`)?.value || '').trim();
      return {
        param_symbol: param,
        param_min: min,
        param_max: max,
        n_points: points,
        output_exprs: expr ? [expr] : []
      };
    }
    case 'rop-cloud-params': {
      const mode = document.getElementById(`${nodeId}-sampling-mode`)?.value || 'x_space';
      const samples = parseInt(document.getElementById(`${nodeId}-samples`)?.value || '10000');
      const span = parseInt(document.getElementById(`${nodeId}-span`)?.value || '6');
      const logxMin = parseFloat(document.getElementById(`${nodeId}-logx-min`)?.value || '-6');
      const logxMax = parseFloat(document.getElementById(`${nodeId}-logx-max`)?.value || '6');
      const targetSpecies = document.getElementById(`${nodeId}-target-species`)?.value || '';
      return { mode, samples, span, logxMin, logxMax, targetSpecies, config: nodeRegistry[nodeId]?.data?.config };
    }
    case 'fret-params': {
      const grid = parseInt(document.getElementById(`${nodeId}-grid`)?.value || '80');
      const min = parseFloat(document.getElementById(`${nodeId}-min`)?.value || '-6');
      const max = parseFloat(document.getElementById(`${nodeId}-max`)?.value || '6');
      return { grid, min, max, config: nodeRegistry[nodeId]?.data?.config };
    }
    case 'scan-2d-params': {
      const param1 = document.getElementById(`${nodeId}-param1`)?.value || '';
      const param2 = document.getElementById(`${nodeId}-param2`)?.value || '';
      const min1 = parseFloat(document.getElementById(`${nodeId}-min1`)?.value || '-6');
      const max1 = parseFloat(document.getElementById(`${nodeId}-max1`)?.value || '6');
      const min2 = parseFloat(document.getElementById(`${nodeId}-min2`)?.value || '-6');
      const max2 = parseFloat(document.getElementById(`${nodeId}-max2`)?.value || '6');
      const points = parseInt(document.getElementById(`${nodeId}-points`)?.value || '50');
      const expr = (document.getElementById(`${nodeId}-expr`)?.value || '').trim();
      return {
        param1_symbol: param1,
        param2_symbol: param2,
        param1_min: min1,
        param1_max: max1,
        param2_min: min2,
        param2_max: max2,
        n_grid: points,
        output_expr: expr
      };
    }
    case 'rop-poly-params': {
      const expr = document.getElementById(`${nodeId}-expr`)?.value || '';
      const param1 = document.getElementById(`${nodeId}-param1`)?.value || '';
      const param2 = document.getElementById(`${nodeId}-param2`)?.value || '';
      const asymptotic = document.getElementById(`${nodeId}-asymptotic`)?.checked ?? true;
      const maxVertices = parseInt(document.getElementById(`${nodeId}-max-vertices`)?.value || '1000');
      return { expr, param1, param2, asymptotic, maxVertices, config: nodeRegistry[nodeId]?.data?.config };
    }
    default:
      return {};
  }
}

function serializeState() {
  const nodes = [];
  for (const [id, info] of Object.entries(nodeRegistry)) {
    const el = document.getElementById(id);
    if (!el) continue;
    nodes.push({
      id,
      type: info.type,
      x: parseFloat(el.style.left) || 0,
      y: parseFloat(el.style.top) || 0,
      width: el.offsetWidth,
      data: getNodeSerialData(id, info.type),
    });
  }
  return {
    version: 1,
    timestamp: new Date().toISOString(),
    canvas: { panX, panY, scale },
    nodes,
    connections: connections.map(c => ({ ...c })),
  };
}

function saveState() {
  const data = serializeState();
  const blob = new Blob([JSON.stringify(data, null, 2)], { type: 'application/json' });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = `rop-workspace-${Date.now()}.json`;
  a.click();
  URL.revokeObjectURL(url);
  showToast('Workspace saved');
}

function loadState() {
  const input = document.createElement('input');
  input.type = 'file';
  input.accept = '.json';
  input.onchange = (e) => {
    const file = e.target.files[0];
    if (!file) return;
    const reader = new FileReader();
    reader.onload = (ev) => {
      try {
        const data = JSON.parse(ev.target.result);
        if (!data.version || !data.nodes) throw new Error('Invalid workspace file');
        applyState(data);
        showToast('Workspace loaded');
      } catch (err) {
        showToast('Failed to load: ' + err.message);
      }
    };
    reader.readAsText(file);
  };
  input.click();
}

function applyState(data) {
  // 1. Clear existing canvas
  for (const id of Object.keys(nodeRegistry)) {
    const el = document.getElementById(id);
    if (el) el.remove();
  }
  Object.keys(nodeRegistry).forEach(k => delete nodeRegistry[k]);
  connections.length = 0;
  svgLayer.querySelectorAll('.wire').forEach(w => w.remove());

  // Reset model state (user will need to rebuild)
  state.sessionId = null;
  state.model = null;
  state.qK_syms = [];

  // 2. Restore canvas position
  if (data.canvas) {
    panX = data.canvas.panX || 0;
    panY = data.canvas.panY || 0;
    scale = data.canvas.scale || 1.0;
    applyViewportTransform();
  }

  // 3. Create nodes and build ID mapping
  const idMap = {}; // oldId → newId
  for (const saved of data.nodes) {
    const newId = createNode(saved.type, saved.x, saved.y);
    if (!newId) continue;
    idMap[saved.id] = newId;

    // Restore width
    const el = document.getElementById(newId);
    if (el && saved.width) el.style.width = `${saved.width}px`;

    // Restore node-specific data
    restoreNodeData(newId, saved.type, saved.data || {});
  }

  // 4. Restore connections (with remapped IDs)
  if (data.connections) {
    for (const conn of data.connections) {
      const fromNode = idMap[conn.fromNode];
      const toNode = idMap[conn.toNode];
      if (fromNode && toNode) {
        connections.push({
          fromNode,
          fromPort: conn.fromPort,
          toNode,
          toPort: conn.toPort,
        });
      }
    }
  }

  // 5. Update wires
  updateConnections();

  // 6. Refresh ROP cloud target options now that connections are restored
  for (const [id, info] of Object.entries(nodeRegistry)) {
    if (info.type !== 'rop-cloud' && info.type !== 'rop-cloud-params') continue;
    updateROPCloudMode(id);
    const savedTarget = info.data?.targetSpecies;
    const sel = document.getElementById(`${id}-target-species`);
    if (sel && savedTarget && Array.from(sel.options).some(o => o.value === savedTarget)) {
      sel.value = savedTarget;
    }
  }
}

function restoreNodeData(nodeId, type, data) {
  switch (type) {
    case 'reaction-network': {
      // Clear default rows added by onInit
      const list = document.getElementById(`${nodeId}-reactions-list`);
      if (list) list.innerHTML = '';
      // Add saved reactions
      if (data.reactions && data.reactions.length > 0) {
        data.reactions.forEach(r => addReactionRow(nodeId, r.rule, r.kd));
      }
      break;
    }
    case 'rop-cloud': {
      const modeEl = document.getElementById(`${nodeId}-sampling-mode`);
      const samplesEl = document.getElementById(`${nodeId}-samples`);
      const spanEl = document.getElementById(`${nodeId}-span`);
      const logxMinEl = document.getElementById(`${nodeId}-logx-min`);
      const logxMaxEl = document.getElementById(`${nodeId}-logx-max`);
      const targetEl = document.getElementById(`${nodeId}-target-species`);
      if (modeEl && data.mode) modeEl.value = data.mode;
      if (samplesEl && data.samples != null) samplesEl.value = data.samples;
      if (spanEl && data.span != null) spanEl.value = data.span;
      if (logxMinEl && data.logxMin != null) logxMinEl.value = data.logxMin;
      if (logxMaxEl && data.logxMax != null) logxMaxEl.value = data.logxMax;
      if (nodeRegistry[nodeId]) nodeRegistry[nodeId].data.targetSpecies = data.targetSpecies || '';
      updateROPCloudMode(nodeId);
      if (targetEl && data.targetSpecies) targetEl.value = data.targetSpecies;
      break;
    }
    case 'fret-heatmap': {
      const gridEl = document.getElementById(`${nodeId}-grid`);
      if (gridEl && data.grid != null) gridEl.value = data.grid;
      break;
    }
    case 'scan-1d-params': {
      const paramEl = document.getElementById(`${nodeId}-param`);
      const minEl = document.getElementById(`${nodeId}-min`);
      const maxEl = document.getElementById(`${nodeId}-max`);
      const pointsEl = document.getElementById(`${nodeId}-points`);
      const exprEl = document.getElementById(`${nodeId}-expr`);
      if (paramEl && data.param_symbol) paramEl.value = data.param_symbol;
      if (minEl && data.param_min != null) minEl.value = data.param_min;
      if (maxEl && data.param_max != null) maxEl.value = data.param_max;
      if (pointsEl && data.n_points != null) pointsEl.value = data.n_points;
      if (exprEl && data.output_exprs && data.output_exprs[0]) exprEl.value = data.output_exprs[0];
      break;
    }
    case 'rop-cloud-params': {
      const modeEl = document.getElementById(`${nodeId}-sampling-mode`);
      const samplesEl = document.getElementById(`${nodeId}-samples`);
      const spanEl = document.getElementById(`${nodeId}-span`);
      const logxMinEl = document.getElementById(`${nodeId}-logx-min`);
      const logxMaxEl = document.getElementById(`${nodeId}-logx-max`);
      const targetEl = document.getElementById(`${nodeId}-target-species`);
      if (modeEl && data.mode) modeEl.value = data.mode;
      if (samplesEl && data.samples != null) samplesEl.value = data.samples;
      if (spanEl && data.span != null) spanEl.value = data.span;
      if (logxMinEl && data.logxMin != null) logxMinEl.value = data.logxMin;
      if (logxMaxEl && data.logxMax != null) logxMaxEl.value = data.logxMax;
      if (nodeRegistry[nodeId]) nodeRegistry[nodeId].data.targetSpecies = data.targetSpecies || '';
      updateROPCloudMode(nodeId);
      if (targetEl && data.targetSpecies) targetEl.value = data.targetSpecies;
      if (data.config && nodeRegistry[nodeId]) {
        nodeRegistry[nodeId].data.config = data.config;
      }
      break;
    }
    case 'fret-params': {
      const gridEl = document.getElementById(`${nodeId}-grid`);
      const minEl = document.getElementById(`${nodeId}-min`);
      const maxEl = document.getElementById(`${nodeId}-max`);
      if (gridEl && data.grid != null) gridEl.value = data.grid;
      if (minEl && data.min != null) minEl.value = data.min;
      if (maxEl && data.max != null) maxEl.value = data.max;
      if (data.config && nodeRegistry[nodeId]) {
        nodeRegistry[nodeId].data.config = data.config;
      }
      break;
    }
    case 'siso-params': {
      const selectEl = document.getElementById(`${nodeId}-siso-select`);
      const minEl = document.getElementById(`${nodeId}-min`);
      const maxEl = document.getElementById(`${nodeId}-max`);
      if (selectEl && data.changeQK) selectEl.value = data.changeQK;
      if (minEl && data.min != null) minEl.value = data.min;
      if (maxEl && data.max != null) maxEl.value = data.max;
      if (data.config && nodeRegistry[nodeId]) {
        nodeRegistry[nodeId].data.config = data.config;
      }
      break;
    }
    case 'scan-2d-params': {
      const param1El = document.getElementById(`${nodeId}-param1`);
      const param2El = document.getElementById(`${nodeId}-param2`);
      const min1El = document.getElementById(`${nodeId}-min1`);
      const max1El = document.getElementById(`${nodeId}-max1`);
      const min2El = document.getElementById(`${nodeId}-min2`);
      const max2El = document.getElementById(`${nodeId}-max2`);
      const pointsEl = document.getElementById(`${nodeId}-points`);
      const exprEl = document.getElementById(`${nodeId}-expr`);
      if (param1El && data.param1_symbol) param1El.value = data.param1_symbol;
      if (param2El && data.param2_symbol) param2El.value = data.param2_symbol;
      if (min1El && data.param1_min != null) min1El.value = data.param1_min;
      if (max1El && data.param1_max != null) max1El.value = data.param1_max;
      if (min2El && data.param2_min != null) min2El.value = data.param2_min;
      if (max2El && data.param2_max != null) max2El.value = data.param2_max;
      if (pointsEl && data.n_grid != null) pointsEl.value = data.n_grid;
      if (exprEl && data.output_expr) exprEl.value = data.output_expr;
      break;
    }
    case 'rop-poly-params': {
      const exprEl = document.getElementById(`${nodeId}-expr`);
      const param1El = document.getElementById(`${nodeId}-param1`);
      const param2El = document.getElementById(`${nodeId}-param2`);
      const asymptoticEl = document.getElementById(`${nodeId}-asymptotic`);
      const maxVerticesEl = document.getElementById(`${nodeId}-max-vertices`);
      if (exprEl && data.expr) exprEl.value = data.expr;
      if (param1El && data.param1) param1El.value = data.param1;
      if (param2El && data.param2) param2El.value = data.param2;
      if (asymptoticEl && data.asymptotic != null) asymptoticEl.checked = data.asymptotic;
      if (maxVerticesEl && data.maxVertices != null) maxVerticesEl.value = data.maxVertices;
      if (data.config && nodeRegistry[nodeId]) {
        nodeRegistry[nodeId].data.config = data.config;
      }
      break;
    }
  }
}

// ===== Parameter Scan 1D Helper Functions =====
const plotResizeObservers = new Map(); // nodeId -> ResizeObserver

function setupPlotResize(nodeId, plotId) {
  // Clean up existing observer
  if (plotResizeObservers.has(nodeId)) {
    plotResizeObservers.get(nodeId).disconnect();
  }

  const plotEl = document.getElementById(plotId);
  if (!plotEl) return;

  const observer = new ResizeObserver(() => {
    Plotly.Plots.resize(plotEl);
  });

  observer.observe(plotEl);
  plotResizeObservers.set(nodeId, observer);
}

function cleanupPlotResize(nodeId) {
  if (plotResizeObservers.has(nodeId)) {
    plotResizeObservers.get(nodeId).disconnect();
    plotResizeObservers.delete(nodeId);
  }
}

function insertSpecies1D(nodeId) {
  const helper = document.getElementById(`${nodeId}-species-helper`);
  const expr = document.getElementById(`${nodeId}-expr`);
  if (helper.value && expr) {
    expr.value += (expr.value ? ' + ' : '') + helper.value;
    helper.value = '';
    const info = nodeRegistry[nodeId];
    if (info && info.type === 'scan-1d-params') {
      triggerConfigUpdate(nodeId, info.type);
    }
  }
}

function updateScan1DConfig(nodeId) {
  const paramSelect = document.getElementById(`${nodeId}-param`);
  const minInput = document.getElementById(`${nodeId}-min`);
  const maxInput = document.getElementById(`${nodeId}-max`);
  const pointsInput = document.getElementById(`${nodeId}-points`);
  const exprInput = document.getElementById(`${nodeId}-expr`);

  if (!paramSelect.value) {
    alert('Please select a scan parameter');
    return;
  }

  if (!exprInput.value.trim()) {
    alert('Please enter an output expression');
    return;
  }

  // Store config in node data
  nodeRegistry[nodeId].data.config = {
    param_symbol: paramSelect.value,
    param_min: parseFloat(minInput.value),
    param_max: parseFloat(maxInput.value),
    n_points: parseInt(pointsInput.value),
    output_exprs: [exprInput.value.trim()],
  };

  showToast('Configuration updated');
}

async function runParameterScan1D(nodeId) {
  const paramSymbol = document.getElementById(`${nodeId}-param`).value;
  const min = parseFloat(document.getElementById(`${nodeId}-min`).value);
  const max = parseFloat(document.getElementById(`${nodeId}-max`).value);
  const points = parseInt(document.getElementById(`${nodeId}-points`).value);

  const expr = document.getElementById(`${nodeId}-expr`).value.trim();
  if (!expr) {
    alert('Please enter an output expression');
    return;
  }

  setNodeLoading(nodeId, true);
  const contentEl = document.getElementById(`${nodeId}-content`);

  try {
    const data = await api('parameter_scan_1d', {
      session_id: state.sessionId,
      param_symbol: paramSymbol,
      param_min: min,
      param_max: max,
      n_points: points,
      output_exprs: [expr],
    });

    contentEl.innerHTML = `<div class="plot-container" id="${nodeId}-plot"></div>`;
    setTimeout(() => {
      plotParameterScan1D(data, `${nodeId}-plot`);
      setupPlotResize(nodeId, `${nodeId}-plot`);
    }, 50);
  } catch (e) {
    contentEl.innerHTML = `<div class="node-error">${e.message}</div>`;
  } finally {
    setNodeLoading(nodeId, false);
  }
}

async function executeScan1DResult(nodeId) {
  // Find connected params node
  const conn = connections.find(c => c.toNode === nodeId && c.toPort === 'params');
  if (!conn) {
    alert('Please connect to a Scan 1D Config node');
    return;
  }

  const paramsNode = nodeRegistry[conn.fromNode];
  if (!paramsNode) {
    alert('Config node has no configuration. Please configure it first.');
    return;
  }

  // Sync latest DOM values to avoid stale config when user changed fields just before Run.
  triggerConfigUpdate(conn.fromNode, paramsNode.type || 'scan-1d-params');
  const config = paramsNode.data?.config;
  if (!config) {
    alert('Config node has no configuration. Please configure it first.');
    return;
  }

  // Validate required fields
  if (!config.param_symbol) {
    alert('Please select a scan parameter in the config node');
    return;
  }

  if (!config.output_exprs || !config.output_exprs[0]) {
    alert('Please enter an output expression in the config node');
    return;
  }

  setNodeLoading(nodeId, true);
  const contentEl = document.getElementById(`${nodeId}-content`);

  try {
    const data = await api('parameter_scan_1d', {
      session_id: state.sessionId,
      param_symbol: config.param_symbol,
      param_min: config.param_min,
      param_max: config.param_max,
      n_points: config.n_points,
      output_exprs: config.output_exprs,
    });

    contentEl.innerHTML = `<div class="plot-container" id="${nodeId}-plot"></div>`;
    setTimeout(() => {
      plotParameterScan1D(data, `${nodeId}-plot`);
      setupPlotResize(nodeId, `${nodeId}-plot`);
    }, 50);
  } catch (e) {
    contentEl.innerHTML = `<div class="node-error">${e.message}</div>`;
  } finally {
    setNodeLoading(nodeId, false);
  }
}

function updateROPCloudConfig(nodeId) {
  const mode = document.getElementById(`${nodeId}-sampling-mode`).value;
  const nSamples = parseInt(document.getElementById(`${nodeId}-samples`).value);

  const config = {
    mode: mode,
    n_samples: nSamples,
  };

  if (mode === 'x_space') {
    const targetSpecies = document.getElementById(`${nodeId}-target-species`)?.value || '';
    const logxMin = parseFloat(document.getElementById(`${nodeId}-logx-min`).value);
    const logxMax = parseFloat(document.getElementById(`${nodeId}-logx-max`).value);
    config.target_species = targetSpecies;
    config.logx_min = logxMin;
    config.logx_max = logxMax;
  } else {
    const span = parseInt(document.getElementById(`${nodeId}-span`).value);
    config.span = span;
  }

  nodeRegistry[nodeId].data.config = config;
  showToast('Configuration updated');
}

async function executeROPCloudResult(nodeId) {
  const conn = connections.find(c => c.toNode === nodeId && c.toPort === 'params');
  if (!conn) {
    alert('Please connect to a ROP Cloud Config node');
    return;
  }

  const paramsNode = nodeRegistry[conn.fromNode];
  if (!paramsNode || !paramsNode.data.config) {
    alert('Config node has no configuration. Please configure it first.');
    return;
  }

  const config = paramsNode.data.config;
  setNodeLoading(nodeId, true);
  const contentEl = document.getElementById(`${nodeId}-content`);

  try {
    let data;
    if (config.mode === 'qk') {
      if (!state.sessionId) throw new Error('Build a model first, or switch to x-space mode');
      const modelConn = connections.find(c => c.toNode === conn.fromNode && c.toPort === 'model');
      if (!modelConn) throw new Error('qK mode requires Model input connection');
      data = await api('rop_cloud', {
        sampling_mode: 'qk',
        session_id: state.sessionId,
        n_samples: config.n_samples,
        span: config.span,
      });
    } else {
      const rxConn = connections.find(c => c.toNode === conn.fromNode && c.toPort === 'reactions');
      if (!rxConn) throw new Error('x-space mode requires Reactions input connection');
      const { reactions } = getReactionsFromNode(rxConn.fromNode);
      if (!reactions.length) throw new Error('Add at least one reaction in the connected Reaction Network');
      data = await api('rop_cloud', {
        sampling_mode: 'x_space',
        reactions: reactions,
        n_samples: config.n_samples,
        logx_min: config.logx_min,
        logx_max: config.logx_max,
        target_species: config.target_species || '',
      });
    }

    contentEl.innerHTML = `<div class="plot-container" id="${nodeId}-plot"></div>`;
    setTimeout(() => {
      plotROPCloud(data, `${nodeId}-plot`);
      setupPlotResize(nodeId, `${nodeId}-plot`);
    }, 50);
  } catch (e) {
    contentEl.innerHTML = `<div class="node-error">${e.message}</div>`;
  } finally {
    setNodeLoading(nodeId, false);
  }
}

function updateFRETConfig(nodeId) {
  const grid = parseInt(document.getElementById(`${nodeId}-grid`).value);
  nodeRegistry[nodeId].data.config = { n_grid: grid };
  showToast('Configuration updated');
}

async function executeFRETResult(nodeId) {
  const conn = connections.find(c => c.toNode === nodeId && c.toPort === 'params');
  if (!conn) {
    alert('Please connect to a FRET Config node');
    return;
  }

  const paramsNode = nodeRegistry[conn.fromNode];
  if (!paramsNode || !paramsNode.data.config) {
    alert('Config node has no configuration. Please configure it first.');
    return;
  }

  const config = paramsNode.data.config;
  setNodeLoading(nodeId, true);
  const contentEl = document.getElementById(`${nodeId}-content`);

  try {
    const data = await api('fret_heatmap', {
      session_id: state.sessionId,
      n_grid: config.n_grid,
    });

    contentEl.innerHTML = `<div class="plot-container" id="${nodeId}-plot"></div>`;
    setTimeout(() => {
      plotHeatmap(data, `${nodeId}-plot`);
      setupPlotResize(nodeId, `${nodeId}-plot`);
    }, 50);
  } catch (e) {
    contentEl.innerHTML = `<div class="node-error">${e.message}</div>`;
  } finally {
    setNodeLoading(nodeId, false);
  }
}

function insertSpecies2D(nodeId) {
  const helper = document.getElementById(`${nodeId}-species-helper`);
  const expr = document.getElementById(`${nodeId}-expr`);
  if (helper.value && expr) {
    expr.value += (expr.value ? ' + ' : '') + helper.value;
    helper.value = '';
    const info = nodeRegistry[nodeId];
    if (info && info.type === 'scan-2d-params') {
      triggerConfigUpdate(nodeId, info.type);
    }
  }
}

function updateScan2DConfig(nodeId) {
  const param1 = document.getElementById(`${nodeId}-param1`).value;
  const param2 = document.getElementById(`${nodeId}-param2`).value;
  const min1 = parseFloat(document.getElementById(`${nodeId}-min1`).value);
  const max1 = parseFloat(document.getElementById(`${nodeId}-max1`).value);
  const min2 = parseFloat(document.getElementById(`${nodeId}-min2`).value);
  const max2 = parseFloat(document.getElementById(`${nodeId}-max2`).value);
  const points = parseInt(document.getElementById(`${nodeId}-points`).value);
  const expr = document.getElementById(`${nodeId}-expr`).value.trim();

  if (!param1 || !param2) {
    alert('Please select both parameters');
    return;
  }

  if (!expr) {
    alert('Please enter an output expression');
    return;
  }

  nodeRegistry[nodeId].data.config = {
    param_symbol_1: param1,
    param_symbol_2: param2,
    param_min_1: min1,
    param_max_1: max1,
    param_min_2: min2,
    param_max_2: max2,
    n_points: points,
    output_exprs: [expr],
  };

  showToast('Configuration updated');
}

async function executeScan2DResult(nodeId) {
  const conn = connections.find(c => c.toNode === nodeId && c.toPort === 'params');
  if (!conn) {
    alert('Please connect to a Scan 2D Config node');
    return;
  }

  const paramsNode = nodeRegistry[conn.fromNode];
  if (!paramsNode) {
    alert('Config node has no configuration. Please configure it first.');
    return;
  }

  // Sync latest DOM values to avoid stale config when user changed fields just before Run.
  triggerConfigUpdate(conn.fromNode, paramsNode.type || 'scan-2d-params');
  const config = paramsNode.data?.config;
  if (!config) {
    alert('Config node has no configuration. Please configure it first.');
    return;
  }

  // Validate required fields
  if (!config.param1_symbol || !config.param2_symbol) {
    alert('Please select both X and Y axis parameters in the config node');
    return;
  }

  if (!config.output_expr) {
    alert('Please enter an output expression in the config node');
    return;
  }

  setNodeLoading(nodeId, true);
  const contentEl = document.getElementById(`${nodeId}-content`);

  try {
    const data = await api('parameter_scan_2d', {
      session_id: state.sessionId,
      ...config
    });

    contentEl.innerHTML = `<div class="plot-container" id="${nodeId}-plot"></div>`;
    setTimeout(() => {
      plotParameterScan2D(data, `${nodeId}-plot`);
      setupPlotResize(nodeId, `${nodeId}-plot`);
    }, 50);
  } catch (e) {
    contentEl.innerHTML = `<div class="node-error">${e.message}</div>`;
  } finally {
    setNodeLoading(nodeId, false);
  }
}

function updateROPPolyConfig(nodeId) {
  const expr = document.getElementById(`${nodeId}-expr`).value.trim();
  const param1 = document.getElementById(`${nodeId}-param1`).value;
  const param2 = document.getElementById(`${nodeId}-param2`).value;
  const asymptotic = document.getElementById(`${nodeId}-asymptotic`).checked;
  const maxVertices = parseInt(document.getElementById(`${nodeId}-max-vertices`).value);

  if (!expr) {
    alert('Please enter an output expression');
    return;
  }

  if (!param1 || !param2) {
    alert('Please select both parameters');
    return;
  }

  nodeRegistry[nodeId].data.config = {
    output_expr: expr,
    param_symbol_1: param1,
    param_symbol_2: param2,
    asymptotic_only: asymptotic,
    max_vertices: maxVertices,
  };

  showToast('Configuration updated');
}

async function executeROPPolyResult(nodeId) {
  const conn = connections.find(c => c.toNode === nodeId && c.toPort === 'params');
  if (!conn) {
    alert('Please connect to a ROP Polyhedron Config node');
    return;
  }

  const paramsNode = nodeRegistry[conn.fromNode];
  if (!paramsNode || !paramsNode.data.config) {
    alert('Config node has no configuration. Please configure it first.');
    return;
  }

  const config = paramsNode.data.config;
  setNodeLoading(nodeId, true);
  const contentEl = document.getElementById(`${nodeId}-content`);

  try {
    const data = await api('rop_polyhedron', {
      session_id: state.sessionId,
      ...config
    });

    contentEl.innerHTML = `<div class="plot-container" id="${nodeId}-plot"></div>`;
    setTimeout(() => {
      plotROPPolyhedron(data, `${nodeId}-plot`);
      setupPlotResize(nodeId, `${nodeId}-plot`);
    }, 50);
  } catch (e) {
    contentEl.innerHTML = `<div class="node-error">${e.message}</div>`;
  } finally {
    setNodeLoading(nodeId, false);
  }
}

function plotParameterScan1D(data, plotId) {
  const { param_symbol, param_values, output_exprs, output_traj } = data;

  const traces = output_exprs.map((expr, i) => ({
    x: param_values,
    y: output_traj.map(row => row[i]),
    mode: 'lines',
    name: expr,
    line: { width: 2 },
  }));

  const layout = {
    autosize: true,
    paper_bgcolor: '#111',
    plot_bgcolor: '#1a1a1a',
    font: { color: '#888', size: 10 },
    margin: { t: 40, b: 60, l: 70, r: 20 },
    title: { text: `Parameter Scan: ${param_symbol}`, font: { color: '#888', size: 11 }, y: 0.98, yanchor: 'top' },
    xaxis: { title: `log10(${param_symbol})`, gridcolor: '#333' },
    yaxis: { title: 'log10(concentration)', gridcolor: '#333' },
    legend: { x: 1, xanchor: 'right', y: 1 },
  };

  Plotly.newPlot(plotId, traces, layout, { responsive: true, displayModeBar: false });
}

// ===== Parameter Scan 2D Helper Functions =====
function insertSpecies2D(nodeId) {
  const helper = document.getElementById(`${nodeId}-species-helper`);
  const expr = document.getElementById(`${nodeId}-expr`);
  if (helper.value && expr) {
    expr.value += (expr.value ? ' + ' : '') + helper.value;
    helper.value = '';
    const info = nodeRegistry[nodeId];
    if (info && info.type === 'scan-2d-params') {
      triggerConfigUpdate(nodeId, info.type);
    }
  }
}

async function runParameterScan2D(nodeId) {
  const param1 = document.getElementById(`${nodeId}-param1`).value;
  const param2 = document.getElementById(`${nodeId}-param2`).value;
  const min1 = parseFloat(document.getElementById(`${nodeId}-min1`).value);
  const max1 = parseFloat(document.getElementById(`${nodeId}-max1`).value);
  const min2 = parseFloat(document.getElementById(`${nodeId}-min2`).value);
  const max2 = parseFloat(document.getElementById(`${nodeId}-max2`).value);
  const grid = parseInt(document.getElementById(`${nodeId}-grid`).value);
  const expr = document.getElementById(`${nodeId}-expr`).value.trim();

  if (!expr) {
    alert('Please enter an output expression');
    return;
  }

  setNodeLoading(nodeId, true);
  const contentEl = document.getElementById(`${nodeId}-content`);

  try {
    const data = await api('parameter_scan_2d', {
      session_id: state.sessionId,
      param1_symbol: param1,
      param2_symbol: param2,
      param1_min: min1,
      param1_max: max1,
      param2_min: min2,
      param2_max: max2,
      n_grid: grid,
      output_expr: expr,
    });

    contentEl.innerHTML = `<div class="plot-container" id="${nodeId}-plot"></div>`;
    setTimeout(() => {
      plotParameterScan2D(data, `${nodeId}-plot`);
      setupPlotResize(nodeId, `${nodeId}-plot`);
    }, 50);
  } catch (e) {
    contentEl.innerHTML = `<div class="node-error">${e.message}</div>`;
  } finally {
    setNodeLoading(nodeId, false);
  }
}

function plotParameterScan2D(data, plotId) {
  const { param1_symbol, param2_symbol, param1_values, param2_values, output_expr, output_grid, regime_grid } = data;

  // Create 3D surface plot
  const traces = [
    {
      z: output_grid,
      x: param1_values,
      y: param2_values,
      type: 'surface',
      colorscale: 'Viridis',
      colorbar: {
        title: `log(${output_expr})`,
     titlefont: { color: '#888', size: 9 },
        tickfont: { size: 8 }
    },
      contours: {
        z: {
        show: true,
          usecolormap: true,
       highlightcolor: "#42f462",
          project: { z: true }
        }
      }
    }
  ];

  const layout = {
    autosize: true,
    paper_bgcolor: '#111',
  plot_bgcolor: '#1a1a1a',
    font: { color: '#888', size: 10 },
    margin: { t: 40, b: 60, l: 70, r: 20 },
    title: {
   text: `${output_expr} vs ${param1_symbol}, ${param2_symbol}`,
      font: { color: '#888', size: 11 },
      y: 0.98,
      yanchor: 'top'
    },
    scene: {
      xaxis: {
        title: `log10(${param1_symbol})`,
        gridcolor: '#333',
        backgroundcolor: '#1a1a1a'
      },
      yaxis: {
        title: `log10(${param2_symbol})`,
        gridcolor: '#333',
        backgroundcolor: '#1a1a1a'
      },
      zaxis: {
        title: `log10(${output_expr})`,
        gridcolor: '#333',
        backgroundcolor: '#1a1a1a'
      },
      bgcolor: '#1a1a1a'
    }
  };

  Plotly.newPlot(plotId, traces, layout, { responsive: true, displayModeBar: false });
}

// ===== ROP Polyhedron Helper Functions =====
function insertSpeciesPoly(nodeId) {
  const helper = document.getElementById(`${nodeId}-species-helper`);
  const expr = document.getElementById(`${nodeId}-expr`);
  if (helper.value && expr) {
    expr.value += (expr.value ? ' + ' : '') + helper.value;
    helper.value = '';
  }
}

async function runROPPolyhedron(nodeId) {
  const expr = document.getElementById(`${nodeId}-expr`).value.trim();
  const param1 = document.getElementById(`${nodeId}-param1`).value;
  const param2 = document.getElementById(`${nodeId}-param2`).value;
  const asymptotic = document.getElementById(`${nodeId}-asymptotic`).checked;
  const maxVertices = parseInt(document.getElementById(`${nodeId}-max-vertices`).value);

  if (!expr) {
    alert('Please enter an output expression');
    return;
  }

  setNodeLoading(nodeId, true);
  const contentEl = document.getElementById(`${nodeId}-content`);

  try {
    const data = await api('rop_polyhedron', {
      session_id: state.sessionId,
      output_expr: expr,
      param1_symbol: param1,
      param2_symbol: param2,
      asymptotic_only: asymptotic,
      max_vertices: maxVertices,
    });

    contentEl.innerHTML = `<div class="plot-container" id="${nodeId}-plot"></div>`;
    setTimeout(() => {
      plotROPPolyhedron(data, `${nodeId}-plot`);
      setupPlotResize(nodeId, `${nodeId}-plot`);
    }, 50);
  } catch (e) {
    contentEl.innerHTML = `<div class="node-error">${e.message}</div>`;
  } finally {
    setNodeLoading(nodeId, false);
  }
}

function plotROPPolyhedron(data, plotId) {
  const { output_expr, param1_symbol, param2_symbol, vertices, edges } = data;

  const traces = [];

  // Add edges
  edges.forEach(edge => {
    traces.push({
      x: edge.ro1,
      y: edge.ro2,
      mode: 'lines',
      line: { color: '#00ff00', width: 2 },
      showlegend: false,
      hoverinfo: 'skip',
    });
  });

  // Add vertices
  if (vertices.length > 0) {
    const vertexRO1 = vertices.map(v => v.ro1);
    const vertexRO2 = vertices.map(v => v.ro2);
    const hoverText = vertices.map(v => `Vertex ${v.idx}<br>Nullity: ${v.nullity}<br>Perm: [${v.perm.join(',')}]`);

    traces.push({
      x: vertexRO1,
      y: vertexRO2,
      mode: 'markers',
      marker: { color: '#ff0000', size: 8 },
      name: 'Vertices',
      hovertext: hoverText,
      hoverinfo: 'text',
    });
  }

  const layout = {
    autosize: true,
    paper_bgcolor: '#111',
    plot_bgcolor: '#1a1a1a',
    font: { color: '#888', size: 10 },
    margin: { t: 40, b: 60, l: 70, r: 20 },
    title: { text: `ROP Polyhedron: ${output_expr}`, font: { color: '#888', size: 11 }, y: 0.98, yanchor: 'top' },
    xaxis: { title: `∂log(${output_expr})/∂log(${param1_symbol})`, gridcolor: '#333' },
    yaxis: { title: `∂log(${output_expr})/∂log(${param2_symbol})`, gridcolor: '#333' },
    showlegend: true,
  };

  Plotly.newPlot(plotId, traces, layout, { responsive: true, displayModeBar: false });
}

// ===== Auto-update Config Nodes =====
function setupAutoUpdate(nodeId, nodeType) {
  const node = document.getElementById(nodeId);
  if (!node) return;

  // Save initial default config
  triggerConfigUpdate(nodeId, nodeType);

  // Find all inputs with auto-update class
  const inputs = node.querySelectorAll('.auto-update');

  inputs.forEach(input => {
    const eventType = input.tagName === 'SELECT' ? 'change' :
                      input.type === 'checkbox' ? 'change' : 'input';

    input.addEventListener(eventType, () => {
      // Debounce for text inputs
      if (input.type === 'text' || input.type === 'number') {
        clearTimeout(input._autoUpdateTimer);
        input._autoUpdateTimer = setTimeout(() => {
          triggerConfigUpdate(nodeId, nodeType);
        }, 500);
      } else {
        // Immediate for selects and checkboxes
        triggerConfigUpdate(nodeId, nodeType);
      }
    });
  });
}

function triggerConfigUpdate(nodeId, nodeType) {
  // Store config in node data
  const info = nodeRegistry[nodeId];
  if (!info) return;

  info.data = info.data || {};
  info.data.config = getNodeSerialData(nodeId, nodeType);
}


// ===== Auto-build Model =====
function setupAutoModelBuild(nodeId) {
  // Check if there's a connection to reaction-network
  const checkAndBuild = () => {
    const rxConn = connections.find(c => c.toNode === nodeId && c.toPort === 'reactions');
    if (rxConn) {
      const { reactions } = getReactionsFromNode(rxConn.fromNode);
      if (reactions.length > 0) {
        // Valid reactions exist, auto-build
        setTimeout(() => buildModel(nodeId), 100);
      }
    }
  };

  // Initial check
  checkAndBuild();

  // Store the check function for later use
  if (!nodeRegistry[nodeId]) return;
  nodeRegistry[nodeId]._autoBuildCheck = checkAndBuild;
}

// Trigger auto-build when reactions change
function triggerAutoModelBuild(reactionNodeId) {
  // Find all connected model-builder nodes
  const modelBuilders = connections
    .filter(c => c.fromNode === reactionNodeId && c.fromPort === 'reactions')
    .map(c => c.toNode);

  modelBuilders.forEach(mbId => {
    const info = nodeRegistry[mbId];
    if (info && info._autoBuildCheck) {
      // Debounce the build
      clearTimeout(info._autoBuildTimer);
      info._autoBuildTimer = setTimeout(() => {
        info._autoBuildCheck();
      }, 500);
    }
  });
}
