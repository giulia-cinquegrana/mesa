Release checklist
=================

This is a guide to what needs to be done before a release can be made.

General steps
-------------

- Pick someone to be a release manager
- Pick a release date 
- Pick a RC1 (release candidate) date which should be ~1 month before the release

Prior to generating a release
-----------------------------

- Update the ZAMS model file by running the work directory found in ``data/star_data/zams_models/create_z2m2_y28``. This may take up to an hour or so. This will generate the file ``data/star_data/zams_models/zams_z2m2_y28.data``. Use the ZAMS model plotting script to verify that the HR diagram and central compositions look reasonable, and commit the new data file.
  
Making a release
----------------

Run the release script in ``MESA_DIR``. This requires ``$MESA_DIR`` to be set and takes one argument the release version (you should add the ``r`` prefix as well).
For version ``r12345`` this script will make a branch ``release/r12345`` and then it:

- Updates :file:`data/version_number`

.. note::
    :file:`data/version_number` is normally not included in a commit (it must be explicitly added via ``git add -f``. as we gitignore it).

- Updates :file:`docs/source/conf.py`
- Updates :file:`Doxyfile`

To the new version ``r12345``.

this script will also make zip archive which can be used for local testing to make sure the release builds.

.. note::
    This zip folder is not what we release. The actual zip folder is generated by Github, so that should be tested as well once it has been made.

The release script does not push any changes to Github. That must be done manually with a git push.


Removing files
--------------

Any files that should not be part of the final release should be added to the ``.gitattributes`` file.
This will prevent the file(s) or folders from appearing in the zip archive.


Documentation
-------------

- The Changelog should be updated.

.. note::
    At a minimum this should mention options that are removed/replaced and how to convert from a previous version to the newest version.

- A release notes document should be written.

- The release branch or tag should be added to the `list of active versions on ReadTheDocs <https://readthedocs.org/projects/mesa-doc/versions/>`__.


Testing
-------

- TestHub should report all tests pass for both Linux and macOS on multiple machines and with different OS versions
- The previous SDK version should be tested.

.. note::
    If the previous SDK does not pass we can decide whether to bump the minimum SDK version or fix the issues.

- A non-SDK machine should test the test_suite.
- At least one Windows machine should get tested.
- Recalibrate test suite cases (things like simplex_solar_calibration and example_astero).


Additional Testing
------------------

Additional checks that are not essential but should be done if there is time.

- Check test_memory runs and reports no memory leaks.
- Check ``MESA`` compiles with ``SHARED_LIBS=True``
- Run with FPE checking on.


Linters
-------

There are a number of linters in the ``linters`` folder. The following MUST be run before release:

- check_photos.py This makes sure the photo read/writes are in sync and is needed to ensure photos work.

.. note::
    If any thing was added or removed from a photo remember to bump the version ``star_def_version`` in ``star_data/public/star_data_def.inc``

- fix_inlists.py This makes sure certain options are disabled in the test suite.

Other linters should be run if possible.


Release steps
-------------

To make an actual release (once testing is complete), first push the git tag made by the release script:

- git push origin release/r12345

This is the key bit, as the Github release will be anchored to this tag.

Goto https://github.com/MESAHub/mesa/releases and craft a new release following the guidelines `here <https://docs.github.com/en/repositories/releasing-projects-on-github/managing-releases-in-a-repository>`_.

.. note::
    If this is a RC release then make sure to click ``This is a pre-release``

Add an appropriate title and description. 

.. note::
    The title should be kept simple like ``Release: r12345``

Once created this zip folder should be downloaded and checked that it installs and runs a test case.

Zenodo
------

Once the zip folder has been created it should be uploaded to Zenodo prior to sending a release announcement.
This helps avoid swamping our GitHub bandwidth with user downloads.

- For a pre-release, do not upload to the main MESA zenodo repository.
  Instead upload to its own Zenodo entry. This can be done on a personal account.
- Official releases need to be uploaded to `this MESA Zenodo page <https://doi.org/10.5281/zenodo.2602941>`_.

Send an email to mesa-users
---------------------------

Send an email announcing the release, this should include:

- Link to Zenodo for dowload (not GitHub).
- A brief summary of the changes
- A link to the Changelog
- Highlight any very disruptive changes that might have occurred
- Any new mesa-developers
- Acknowledge those in the community who have helped in some way during this release (bug reports, PR's, testing during the RC phase, being very active on mesa-users)
- Remind people that we welcome any contributions (big or small)

Acknowledging support
---------------------

Getting all authors who committed code (this includes merged pull requests) ::

    git log --format='%aN' r21.12.1..HEAD | sort -u


Listing all commits that acknowledge help from someone ::

    git log --all --grep="-by" r21.12.1..HEAD



Post release fixes
------------------

By having the release be in a separate branch we can push changes if we need to to fix issues however this should be done with caution. Changes to the documentation (highlighting some workaround
are fine). Making changes to the code itself is more tricky (due to the Zenodo upload being fixed and change requiring a new Zenodo upload). It may be easier if a version
needs fixes to simply push a new release, and flag the current release as not working.

New readthedocs version
-----------------------
 
First gain access to the readthedocs account (that is currently accessible by Rich and Rob). Then:

- Goto the ``Versions`` page
- Find the release branch (not the tag) and ``Activate`` it
- We want the branch not the tag so that we can update the docs post release.
- Wait for it to build and check it works
- Goto ``Admin`` page and then the ``Advanced settings`` tab
- Switch the default version to the release.
- Click ``save`` at the bottom of the page
