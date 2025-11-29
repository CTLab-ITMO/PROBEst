document.addEventListener('DOMContentLoaded', function() {
    const form = document.getElementById('probeForm');
    const submitBtn = document.getElementById('submitBtn');
    const loadingOverlay = document.getElementById('loadingOverlay');
    const resultsSection = document.getElementById('resultsSection');
    const resultsContent = document.getElementById('resultsContent');
    const downloadBtn = document.getElementById('downloadBtn');

    form.addEventListener('submit', async function(e) {
        e.preventDefault();
        
        // Show loading overlay
        loadingOverlay.style.display = 'flex';
        resultsSection.style.display = 'none';
        submitBtn.disabled = true;

        // Create FormData
        const formData = new FormData(form);

        try {
            const response = await fetch('/process', {
                method: 'POST',
                body: formData
            });

            const data = await response.json();

            if (!response.ok) {
                throw new Error(data.error || 'An error occurred');
            }

            if (data.success) {
                displayResults(data.top_probes, data.total_probes);
            } else {
                throw new Error(data.error || 'Processing failed');
            }
        } catch (error) {
            showError(error.message);
        } finally {
            loadingOverlay.style.display = 'none';
            submitBtn.disabled = false;
        }
    });

    downloadBtn.addEventListener('click', async function() {
        try {
            const response = await fetch('/download');
            
            if (!response.ok) {
                const error = await response.json();
                throw new Error(error.error || 'Download failed');
            }

            // Get filename from Content-Disposition header or use default
            const contentDisposition = response.headers.get('Content-Disposition');
            let filename = 'PROBEst_results.zip';
            if (contentDisposition) {
                const filenameMatch = contentDisposition.match(/filename="?(.+)"?/);
                if (filenameMatch) {
                    filename = filenameMatch[1];
                }
            }

            // Create blob and download
            const blob = await response.blob();
            const url = window.URL.createObjectURL(blob);
            const a = document.createElement('a');
            a.href = url;
            a.download = filename;
            document.body.appendChild(a);
            a.click();
            window.URL.revokeObjectURL(url);
            document.body.removeChild(a);
        } catch (error) {
            alert('Error downloading results: ' + error.message);
        }
    });

    function displayResults(topProbes, totalProbes) {
        resultsContent.innerHTML = '';

        // Add summary
        const summary = document.createElement('div');
        summary.className = 'success-message';
        summary.innerHTML = `<strong>Success!</strong> Generated ${totalProbes} probes. Top 5 probes are shown below.`;
        resultsContent.appendChild(summary);

        // Display top probes
        if (topProbes && topProbes.length > 0) {
            topProbes.forEach((probe, index) => {
                const card = createProbeCard(probe, index + 1);
                resultsContent.appendChild(card);
            });
        } else {
            const noProbes = document.createElement('div');
            noProbes.className = 'error-message';
            noProbes.textContent = 'No probes found in results.';
            resultsContent.appendChild(noProbes);
        }

        // Show results section
        resultsSection.style.display = 'block';
        resultsSection.scrollIntoView({ behavior: 'smooth' });
    }

    function createProbeCard(probe, rank) {
        const card = document.createElement('div');
        card.className = 'probe-card';

        const title = document.createElement('h3');
        title.textContent = `#${rank} - ${probe.name || 'Unnamed Probe'}`;
        card.appendChild(title);

        const info = document.createElement('div');
        info.className = 'probe-info';

        const hitsItem = document.createElement('div');
        hitsItem.className = 'info-item';
        hitsItem.innerHTML = `
            <span class="info-label">Hits</span>
            <span class="info-value">${probe.hits}</span>
        `;
        info.appendChild(hitsItem);

        const lengthItem = document.createElement('div');
        lengthItem.className = 'info-item';
        lengthItem.innerHTML = `
            <span class="info-label">Length</span>
            <span class="info-value">${probe.sequence ? probe.sequence.length : 0} bp</span>
        `;
        info.appendChild(lengthItem);

        card.appendChild(info);

        const sequenceDiv = document.createElement('div');
        sequenceDiv.className = 'probe-sequence';
        sequenceDiv.textContent = probe.sequence || 'No sequence available';
        card.appendChild(sequenceDiv);

        return card;
    }

    function showError(message) {
        resultsContent.innerHTML = '';
        const errorDiv = document.createElement('div');
        errorDiv.className = 'error-message';
        errorDiv.innerHTML = `<strong>Error:</strong> ${message}`;
        resultsContent.appendChild(errorDiv);
        resultsSection.style.display = 'block';
        resultsSection.scrollIntoView({ behavior: 'smooth' });
    }
});


