## Rstudio - Getting started

Connect to RStudio by setting your password (note, password will not be visible on the screen):

```
sudo passwd $USER
```

figuring out your username:

```
echo My username is $USER
```

and finding YOUR RStudio server interface Web address:

```
echo http://$(hostname):8787/
```

Now go to that Web address in your Web browser, and log in with the username and password from above.
