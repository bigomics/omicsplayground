# INSTALLATION


## Installing Orca

Orca is an Electron app that generates images and reports of Plotly
things like plotly.js graphs, dash apps, dashboards from the command
line. Orca is used to export plotly images to PDF/PNG.

	apt install libgtk2.0-0 libgconf-2-4
	npm install -g electron@1.8.4 orca --unsafe-perm=true --allow-root
		
	

## Installing NGINX reverse proxy



- Create a file /etc/systemd/system/docker.service.d/override.conf
