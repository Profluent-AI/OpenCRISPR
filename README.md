![header](imgs/header.png)

# OpenCRISPR

This repository contains releases for OpenCRISPR, a set of free and open gene editing systems designed by Profluent Bio.

## Releases

| Release | Date | Description |
| :-------------- | :------ | :------- |
| [OpenCRISPR-1](OpenCRISPR-1) | 2024-04-22 | AI-designed, RNA-programmable gene editor with NGG PAM preference. <br>Described in [Ruffolo, Nayfach, Gallagher, and Bhatnagar et al., _Nature_ (2025)]([https://www.biorxiv.org/content/10.1101/2024.04.22.590591v1](https://www.nature.com/articles/s41586-025-09298-z)). |

## FAQs

**What is OpenCRISPR-1?** OpenCRISPR-1 is an AI-created gene editor, consisting of a Cas9-like protein and guide RNA, fully developed using Profluent’s large language models (LLMs). The OpenCRISPR-1 protein maintains the prototypical architecture of a Type II Cas9 nuclease but is hundreds of mutations away from SpCas9 or any other known natural CRISPR-associated protein. You can view OpenCRISPR-1 as a drop-in replacement for many protocols that need a cas9-like protein with an NGG PAM and you can even use it with canonical SpCas9 gRNAs. OpenCRISPR-1 can be fused in a deactivated or nickase format for next generation gene editing techniques like base, prime, or epigenome editing. Find out more in our preprint.

**Why are you releasing OpenCRISPR free of charge – what’s the catch?** There is no catch. OpenCRISPR is free for commercial use to any users who take a license. In a world where gene editing technologies can be difficult to access for both researchers and patients for various reasons, we felt the need to put our company mission into action and release some of the byproducts of our prolific protein design engine to enable more discoveries in the gene editing industry. For partners where further customization and expanded features for OpenCRISPR or another system might be desired, we offer a high-touch collaboration model.

**Are you really not asking for anything?** In addition to abiding by our terms of use, we kindly ask that you allow us to acknowledge you as a user and to let us know when any products using OpenCRISPR advance to the clinic or commercial stages.

**Have you filed IP on OpenCRISPR?** Yes.

**If OpenCRISPR is truly open source, then why do I need to sign a license agreement?** The sequence is freely available via the pre-print. We considered many factors to make accessing OpenCRISPR as frictionless and lightweight as possible; chief among these was ensuring its ethical and safe use. For this reason, if OpenCRISPR users wish to use the molecule for commercial therapeutic uses, we require them to execute a simple license agreement that includes obligations to use the tool for ethical purposes only, in addition to other terms of use.

**What does the license include?** The current release includes the protein sequence of OpenCRISPR-1 along with a compatible AI-generated gRNA, though it is also compatible with canonical Cas9 gRNAs.

**Will there be additional OpenCRISPR releases in the future?** Stay tuned…

**Do you provide protocols?** Please see our pre-print in bioRxiv for a general protocol in addition to a readme protocol that accompanies the sequence release. Other general protocols for editing enzymes should also be compatible.

**Is there a way to share my experience using OpenCRISPR with Profluent?** We expressly welcome any feedback on OpenCRISPR and especially sharing of any observations as you’re using the system. If you find that certain attributes could be changed or improved for your particular needs, please reach out!

**OpenCRISPR is interesting, but I have more needs; what does Profluent offer?** We are open to collaboratively iterate and customize an AI-designed solution that is a perfect match for your specific therapeutic application. This ranges from customized gene editors, antibodies, and broader enzymes. Please email `partnerships@profluent.bio`.

## License

OpenCRISPR is free and public for your research and commercial usage. To ensure the ethical and safe commercial use, we have a simple license agreement that includes obligations to use the tool for ethical purposes only, in addition to other terms. Please complete this [form](https://docs.google.com/forms/d/1h3UbiwBgSUJMgR_6o2WlfEvewfE1Ldmar_FrNyazSv4) to gain access to relevant documents and next steps.

## Citing OpenCRISPR

If you use OpenCRISPR in your research, please cite the following preprint:

```bibtex
@article{ruffolo2025design,
  title={Design of highly functional genome editors by modelling CRISPR--Cas sequences},
  author={Ruffolo, Jeffrey A and Nayfach, Stephen and Gallagher, Joseph and Bhatnagar, Aadyot and Beazer, Joel and Hussain, Riffat and Russ, Jordan and Yip, Jennifer and Hill, Emily and Pacesa, Martin and others},
  journal={Nature},
  pages={1--8},
  year={2025},
  publisher={Nature Publishing Group UK London}
}
```
